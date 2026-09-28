import numpy as np

def fits_hdus(fn):
    raw = open(fn,'rb').read()
    pos = 0
    hdus = []
    while pos + 2880 <= len(raw):
        # gather header blocks
        hdr_bytes = b''
        while True:
            block = raw[pos:pos+2880]
            if len(block) < 2880: break
            hdr_bytes += block
            pos += 2880
            if b'END' in block: break
        s = hdr_bytes.decode('ascii','replace')
        cards = [s[i:i+80] for i in range(0, len(s), 80)]
        hdr = {}
        for card in cards:
            c = card.strip()
            if not c or c.startswith(('END','COMMENT','HISTORY')): continue
            key = c[:8].strip()
            val = c[10:].split('/')[0].strip()
            if val:
                try: hdr[key] = int(val)
                except ValueError:
                    try: hdr[key] = float(val)
                    except ValueError: hdr[key] = val.strip("'")
        hdus.append((hdr, pos))
        # data size
        naxis = hdr.get('NAXIS',0)
        if naxis == 0:
            continue
        if hdr.get('XTENSION','').strip() == 'BINTABLE':
            nrows = hdr.get('NAXIS2',0); rowlen = hdr.get('NAXIS1',0)
            pos += nrows*rowlen
            # pad to 2880
            pos += (2880 - pos % 2880) % 2880
        else:
            axes = [hdr.get('NAXIS%d'%(i+1),1) for i in range(naxis)]
            bp = hdr.get('BITPIX',-32)
            sz = 4 if abs(bp)==32 else 8
            pos += int(np.prod(axes))*sz
            pos += (2880 - pos % 2880) % 2880
    return hdus

def fits_get_image(fn):
    """Return (header, data) of the first NAXIS>0 IMAGE or the BINTABLE columns."""
    raw = open(fn,'rb').read()
    hdus = fits_hdus(fn)
    for hdr, pos in hdus:
        if hdr.get('NAXIS',0) == 0: continue
        if hdr.get('XTENSION','').strip() == 'BINTABLE':
            nrows = hdr.get('NAXIS2',0); rowlen = hdr.get('NAXIS1',0)
            nfields = hdr.get('TFIELDS',0)
            cols = []
            off = pos
            for i in range(1, nfields+1):
                fmt = hdr.get('TFORM%d'%i,'')
                n = hdr.get('TDIM%d'%i,'')
                # Planck: each column is 1-D of length NAXIS2? or 2D; handle common cases
                count = 1
                for t in fmt[1:]:
                    if t.isdigit(): count = count*10 + int(t)
                    else: break
                cols.append(count)
            tot = sum(cols)
            dt = '>f8' if 'D' in hdr.get('TFORM1','') else '>f4'
            data = np.frombuffer(raw, dtype=dt, count=tot, offset=pos)
            return hdr, data
        else:
            axes = [hdr.get('NAXIS%d'%(i+1),1) for i in range(hdr.get('NAXIS',0))]
            bp = hdr.get('BITPIX',-32)
            dt = '>f8' if bp == -64 else '>f4'
            data = np.frombuffer(raw, dtype=dt, count=int(np.prod(axes)), offset=pos)
            return hdr, data
    return None, None

if __name__ == '__main__':
    fn = '/Users/admin/Documents/BlackHoles-Infinity/collected-data/smica_2048.fits'
    hdr, d = fits_get_image(fn)
    print("SMICA keys:", {k:hdr[k] for k in ['NAXIS','NAXIS1','NAXIS2','TFIELDS','TFORM1','TFORM2','TFORM3','NSIDE','ORDERING','FIRSTPIX','LASTPIX'] if k in hdr})
    print("data shape:", d.shape, "dtype:", d.dtype, "min/mean/max:", d.min(), d.mean(), d.max())
    fn2 = '/Users/admin/Documents/BlackHoles-Infinity/collected-data/mask_common.fits'
    hdr2, m = fits_get_image(fn2)
    print("MASK:", {k:hdr2[k] for k in ['NAXIS','NAXIS1','NAXIS2','BITPIX','NSIDE','ORDERING'] if k in hdr2})
    print("mask shape:", m.shape, "unique:", np.unique(m)[:6], "sky frac:", m.mean())
