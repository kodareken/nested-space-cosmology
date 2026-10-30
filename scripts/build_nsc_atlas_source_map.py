"""Index every recorded/cited arXiv work without claiming a new full reading."""
import argparse
from collections import defaultdict
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
START='<!-- nsc-source-catalog:start -->'
END='<!-- nsc-source-catalog:end -->'
RX=re.compile(r'https?://(?:www\.)?arxiv\.org/(?:abs|html|pdf)/([a-z-]+/[0-9]{7}|[0-9]{4}\.[0-9]{4,5})(v[0-9]+)?')


def extract(text,kind):
    rows={};heading='Recorded research';section=''
    for line in text.splitlines():
        if kind=='library' and line.startswith('## '):heading=line[3:]
        if kind=='paper' and line.startswith('### '):
            heading=line[4:];match=re.match(r'([A-Z]?\.?\d+(?:\.\d+)*(?:\.)?)\s',heading)
            section=match.group(1).rstrip('.') if match else heading[:45]
        for match in RX.finditer(line):
            ident=match.group(1)
            row=rows.setdefault(ident,dict(group=heading,version=match.group(2) or '',locations=set(),role=''))
            if kind=='library' and line.startswith('- [') and ' — ' in line and not row['role']:
                role=line.split(' — ',1)[1]
                role=re.sub(r'\[([^\]]+)\]\([^)]+\)',r'\1',role).replace('|','/').replace('`','')
                row['role']=role if len(role)<=300 else role[:300].rsplit(' ',1)[0]+' ...'
            if kind=='paper':row['locations'].add(section or heading)
    return rows


def outputs():
    paper=(ROOT/'paper/nested-space-cosmology.md').read_text(encoding='utf-8')
    if START in paper:
        if END not in paper:raise ValueError('incomplete generated catalog marker')
        base=paper.split(START)[0].rstrip()+'\n'
    else:base=paper.rstrip()+'\n'
    lib=extract((ROOT/'docs/arxiv-studies.md').read_text(encoding='utf-8'),'library')
    main=extract(base,'paper');ids=sorted(set(lib)|set(main));groups=defaultdict(list)
    doc=['# Atlas source connections','',
         f'This navigation index contains **{len(ids)} distinct arXiv works**: {len(main)} are cited in the main argument and {len(set(lib)-set(main))} are retained in the reading library only.',
         '', 'The entries locate existing readings and mathematical use. They do not claim that every paper has been read in full, independently verified, or joined into one compatible action. Versions and reading depth remain in [the library](arxiv-studies.md). The manuscript and scientific owners supply the actual connecting equations.',
         '', '| arXiv work | Recorded reading route | Recorded mathematical role | Main manuscript locations |','| --- | --- | --- | --- |']
    for ident in ids:
        row=lib.get(ident,main.get(ident));group=lib[ident]['group'] if ident in lib else 'Main-manuscript references'
        groups[group].append(ident)
        link=f'[{ident}](https://arxiv.org/abs/{ident}{row["version"]})'
        locations=', '.join(sorted(main.get(ident,{}).get('locations',set()))) or 'Library only'
        role=lib.get(ident,{}).get('role') or 'See the linked reading route and manuscript section for the stated use and assumptions.'
        doc.append(f'| {link} | {group.replace("|","/")} | {role} | {locations} |')
    appendix=[START,'## Appendix C. Complete arXiv source route index','',
        f'The index contains {len(ids)} distinct arXiv works recorded in the research library or cited above. It is a map of sources and reading routes, not a claim that all actions or domains have already been matched.',
        '', 'A dagger marks a work retained in the reading library but not cited in the main argument before this appendix. The [detailed index](../docs/nsc-atlas-source-map.md) gives manuscript locations; the [reading library](../docs/arxiv-studies.md) retains versions, inspected passages and unread boundaries.',
        '', '| Recorded reading route | Source identifiers |','| --- | --- |']
    emitted=[]
    for group,items in groups.items():
        for offset in range(0,len(items),5):
            links=[]
            for ident in items[offset:offset+5]:
                row=lib.get(ident,main.get(ident));label=ident+(' †' if ident not in main else '')
                links.append(f'[{label}](https://arxiv.org/abs/{ident}{row["version"]})');emitted.append(ident)
            appendix.append(f'| {group.replace("|","/")} | '+', '.join(links)+' |')
    assert sorted(emitted)==ids and len(emitted)==len(set(emitted))
    appendix.extend(['',END,''])
    return {'docs/nsc-atlas-source-map.md':'\n'.join(doc)+'\n',
            'paper/nested-space-cosmology.md':base+'\n'+'\n'.join(appendix)},dict(unique=len(ids),main_argument=len(main),library_only=len(set(lib)-set(main)),groups=len(groups))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    content,counts=outputs()
    for name,text in content.items():
        path=ROOT/name
        if args.check:
            if not path.exists() or path.read_text(encoding='utf-8')!=text:raise ValueError('source index stale: '+name)
        else:path.write_text(text,encoding='utf-8',newline='\n')
    print('Source index checked' if args.check else 'Source index written',counts)


if __name__=='__main__':main()
