"""One Dirac driver: a lock, a process pool, and sorted family assembly.

Each worker is one family. Assembly order is the family key, not the
completion order. A second derive_nsc_ Dirac process, or a lock whose
pid is still alive, refuses the start.
"""
from concurrent.futures import ProcessPoolExecutor
import os
from pathlib import Path
import subprocess

DEFAULT_WORKERS = 9


def dirac_process_lines():
    if os.name == 'nt':
        # Enumerate only scientific drivers; POSIX ps is unavailable here.
        command = ("Get-CimInstance Win32_Process | Where-Object { "
                   "$_.CommandLine -match 'python' -and "
                   "$_.CommandLine -match 'derive_nsc_' -and "
                   "$_.CommandLine -match '--run' } | ForEach-Object { "
                   "[string]$_.ProcessId + ' ' + $_.CommandLine }")
        output = subprocess.check_output(
            ['powershell.exe', '-NoProfile', '-NonInteractive', '-Command', command],
            text=True, timeout=15)
        return [line.strip() for line in output.splitlines() if line.strip()]
    output = subprocess.check_output(
        ['ps', '-ax', '-o', 'pid=', '-o', 'command='], text=True)
    return [line.strip() for line in output.splitlines() if line.strip()]


def other_dirac_runs(lines, self_pid):
    """A Dirac run is a python derive_nsc_ process started with --run."""
    parent = os.getppid()
    others = []
    for line in lines:
        pid_text, _, command = line.partition(' ')
        try:
            pid = int(pid_text)
        except ValueError:
            continue
        if pid in (int(self_pid), parent):
            continue
        if 'python' in command and 'derive_nsc_' in command and '--run' in command:
            others.append(line)
    return others


def pid_is_alive(pid):
    if os.name == 'nt':
        # os.kill(pid, 0) sends CTRL_C_EVENT on Windows. Query a process handle.
        import ctypes
        from ctypes import wintypes
        native = ctypes.WinDLL('kernel32', use_last_error=True)
        native.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
        native.OpenProcess.restype = wintypes.HANDLE
        native.WaitForSingleObject.argtypes = (wintypes.HANDLE, wintypes.DWORD)
        native.WaitForSingleObject.restype = wintypes.DWORD
        native.CloseHandle.argtypes = (wintypes.HANDLE,)
        native.CloseHandle.restype = wintypes.BOOL
        number = int(pid)
        if number < 1:
            return False
        handle = native.OpenProcess(0x00100000, False, number)  # SYNCHRONIZE
        if not handle:
            error = ctypes.get_last_error()
            if error == 87:  # nonexistent process
                return False
            if error == 5:  # a protected live process
                return True
            raise ctypes.WinError(error)
        try:
            status = native.WaitForSingleObject(handle, 0)
            if status == 258:  # WAIT_TIMEOUT: process has not exited
                return True
            if status == 0:
                return False
            raise ctypes.WinError(ctypes.get_last_error())
        finally:
            native.CloseHandle(handle)
    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def acquire_dirac_lock(lock_path, owner, self_pid, lines=None):
    """Create the lock or raise. A dead pid does not keep the lock."""
    path = Path(lock_path)
    others = other_dirac_runs(
        dirac_process_lines() if lines is None else lines, self_pid)
    if others:
        raise RuntimeError('another derive_nsc_ Dirac run is live: ' + others[0])
    if path.exists():
        recorded = path.read_text().splitlines()
        holder = recorded[0] if recorded else ''
        if holder and pid_is_alive(holder) and int(holder) != int(self_pid):
            raise RuntimeError('Dirac lock is held by live pid ' + holder)
        path.unlink()
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
    descriptor = os.open(path, flags, 0o644)
    try:
        os.write(descriptor, f'{int(self_pid)}\n{owner}\n'.encode())
    finally:
        os.close(descriptor)
    return path


def release_dirac_lock(lock_path, self_pid):
    path = Path(lock_path)
    if not path.exists():
        return
    recorded = path.read_text().splitlines()
    if recorded and recorded[0] == str(int(self_pid)):
        path.unlink()


def orphan_family_payloads(directory):
    directory = Path(directory)
    if not directory.exists():
        return ()
    json_stems = {path.stem for path in directory.glob('family-*.json')}
    npz_stems = {path.stem for path in directory.glob('family-*.npz')}
    orphans = tuple(sorted(json_stems.symmetric_difference(npz_stems)))
    if orphans:
        raise ValueError('orphan family payload: ' + orphans[0])
    return orphans


def echo_item(item):
    return item


def _worker_entry(function, item):
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['OMP_NUM_THREADS'] = '1'
    return function(item)


def map_families(function, items, workers=DEFAULT_WORKERS):
    """Run items in a pool and return results sorted by the submitted order."""
    items = tuple(items)
    if not items:
        return ()
    count = min(int(workers), len(items))
    if count < 1:
        raise ValueError('at least one family worker required')
    with ProcessPoolExecutor(max_workers=count) as pool:
        results = list(pool.map(_worker_entry, [function] * len(items), items))
    return tuple(results)
