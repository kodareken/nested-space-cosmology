"""Real filesystem publication contracts, including the native Windows path."""

from concurrent.futures import ThreadPoolExecutor
import hashlib
import threading

import pytest

from recursive_horizons.evidence_io import (
    PostpublicationUncertainty, PrepublicationError, PublicationReceipt,
    publish_exclusive_directory, publish_exclusive_file,
)


def test_nested_unicode_file_and_directory_are_published_exactly(tmp_path):
    parent = tmp_path / "evidence"
    parent.mkdir()
    payload = b"native publication\0\xff"
    receipt = publish_exclusive_file(tmp_path, "evidence/result-λ.bin", payload)
    assert receipt == PublicationReceipt("evidence/result-λ.bin", "file", hashlib.sha256(payload).hexdigest())
    assert (parent / "result-λ.bin").read_bytes() == payload
    receipt = publish_exclusive_directory(tmp_path, "evidence/terminal",
                                          {"a.txt": b"A", "nested/b.txt": b"B"})
    assert receipt.kind == "directory"
    assert (parent / "terminal/nested/b.txt").read_bytes() == b"B"
    assert not list(parent.glob(".*.evidence-io-stage-*"))


@pytest.mark.parametrize("existing", [b"", b"original", b"new"])
def test_existing_file_is_never_replaced_even_with_identical_bytes(tmp_path, existing):
    target = tmp_path / "saved.bin"
    target.write_bytes(existing)
    with pytest.raises(PrepublicationError) as error:
        publish_exclusive_file(tmp_path, "saved.bin", b"new")
    assert error.value.published is False
    assert target.read_bytes() == existing
    assert not list(tmp_path.glob(".*.evidence-io-stage-*"))


def test_destination_created_at_publication_boundary_is_retained(tmp_path):
    def concurrent_writer(cut):
        if cut == "before_publish":
            (tmp_path / "saved.bin").write_bytes(b"concurrent writer")

    with pytest.raises(PrepublicationError) as error:
        publish_exclusive_file(tmp_path, "saved.bin", b"ours", _fault_hook=concurrent_writer)
    assert error.value.published is False
    assert (tmp_path / "saved.bin").read_bytes() == b"concurrent writer"
    assert not list(tmp_path.glob(".*.evidence-io-stage-*"))


def test_eight_concurrent_native_publishers_have_exactly_one_winner(tmp_path):
    barrier = threading.Barrier(8, timeout=15)

    def worker(number):
        payload = f"writer {number}".encode()

        def ready(cut):
            if cut == "before_publish":
                barrier.wait()

        try:
            return payload, publish_exclusive_file(tmp_path, "shared.bin", payload, _fault_hook=ready)
        except PrepublicationError as error:
            return payload, error

    with ThreadPoolExecutor(max_workers=8) as pool:
        outcomes = list(pool.map(worker, range(8)))
    winners = [(payload, receipt) for payload, receipt in outcomes if isinstance(receipt, PublicationReceipt)]
    assert len(winners) == 1
    assert sum(isinstance(receipt, PrepublicationError) for _, receipt in outcomes) == 7
    assert (tmp_path / "shared.bin").read_bytes() == winners[0][0]
    assert winners[0][1].payload_sha256 == hashlib.sha256(winners[0][0]).hexdigest()
    assert not list(tmp_path.glob(".*.evidence-io-stage-*"))


def test_staging_tamper_fails_before_publication_and_is_retained(tmp_path):
    def tamper(cut):
        if cut == "before_publish":
            next(tmp_path.glob(".saved.bin.evidence-io-stage-*")).write_bytes(b"tampered")

    with pytest.raises(PrepublicationError) as error:
        publish_exclusive_file(tmp_path, "saved.bin", b"ours", _fault_hook=tamper)
    assert error.value.published is False
    assert error.value.staging_retained is True
    assert not (tmp_path / "saved.bin").exists()
    assert next(tmp_path.glob(".saved.bin.evidence-io-stage-*")).read_bytes() == b"tampered"


def test_postpublication_failure_preserves_the_visible_evidence(tmp_path):
    def stop_after_publish(cut):
        if cut == "after_publish":
            raise OSError("postpublication failure")

    with pytest.raises(PostpublicationUncertainty) as error:
        publish_exclusive_file(tmp_path, "saved.bin", b"ours", _fault_hook=stop_after_publish)
    assert error.value.published is True
    assert (tmp_path / "saved.bin").read_bytes() == b"ours"
    with pytest.raises(PrepublicationError):
        publish_exclusive_file(tmp_path, "saved.bin", b"ours")


def test_existing_directory_is_never_replaced(tmp_path):
    target = tmp_path / "terminal"
    target.mkdir()
    (target / "original").write_bytes(b"original")
    with pytest.raises(PrepublicationError):
        publish_exclusive_directory(tmp_path, "terminal", {"result": b"new"})
    assert (target / "original").read_bytes() == b"original"
    assert not (target / "result").exists()


def test_concurrent_directory_publication_retains_one_complete_inventory(tmp_path):
    barrier = threading.Barrier(4, timeout=15)

    def worker(number):
        payload = f"writer {number}".encode()

        def ready(cut):
            if cut == "before_publish":
                barrier.wait()

        try:
            receipt = publish_exclusive_directory(
                tmp_path, "terminal", {"nested/result": payload}, _fault_hook=ready)
            return payload, receipt
        except PrepublicationError as error:
            return payload, error

    with ThreadPoolExecutor(max_workers=4) as pool:
        outcomes = list(pool.map(worker, range(4)))
    winners = [payload for payload, result in outcomes if isinstance(result, PublicationReceipt)]
    assert len(winners) == 1
    assert (tmp_path / "terminal/nested/result").read_bytes() == winners[0]
    assert sum(isinstance(result, PrepublicationError) for _, result in outcomes) == 3
    assert not list(tmp_path.glob(".*.evidence-io-stage-*"))
