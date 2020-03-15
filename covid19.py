"""
Multiple Sequence Alignment Helper and Visualizer

Author: Rick Kim
Email: rkim@insilico.us.com
"""

from aiohttp import web
import os
import subprocess
import argparse
import sys
import webbrowser
import platform
import asyncio

canonical_protein_names = [
    'orf1a polyprotein',
    'orf1ab polyprotein',
    'leader protein', 
    'nsp2', 
    'nsp3', 
    'nsp4', 
    '3c-like proteinase', 
    'nsp6', 
    'nsp7', 
    'nsp8', 
    'nsp9', 
    'nsp10', 
    'rna-dependent rna polymerase', 
    'helicase', 
    '3\'-to-5\' exonuclease', 
    'endornase', 
    '2\'-o-ribose methyltransferase', 
    'surface glycoprotein', 
    'orf3a protein', 
    'envelope protein', 
    'membrane glycoprotein', 
    'orf6 protein', 
    'orf7a protein', 
    'orf7b', 
    'orf8 protein', 
    'nucleocapsid phosphoprotein', 
    'orf10 protein']

msa_grouping = {
    'leader protein': 'orf1ab polyprotein',
    'nsp2': 'orf1ab polyprotein',
    'nsp3': 'orf1ab polyprotein',
    'nsp4': 'orf1ab polyprotein',
    '3c-like proteinase': 'orf1ab polyprotein',
    'nsp6': 'orf1ab polyprotein',
    'nsp7': 'orf1ab polyprotein',
    'nsp8': 'orf1ab polyprotein',
    'nsp9': 'orf1ab polyprotein',
    'nsp10': 'orf1ab polyprotein',
    'rna-dependent rna polymerase': 'orf1ab polyprotein',
    'helicase': 'orf1ab polyprotein',
    '3\'-to-5\' exonuclease': 'orf1ab polyprotein',
    'endornase': 'orf1ab polyprotein',
    '2\'-o-ribose methyltransferase': 'orf1ab polyprotein',
    'orf1a': 'orf1ab polyprotein',
}

def fasta_parser (f):
    while True:
        pos = f.tell()
        line = f.readline()
        if line[0] == '>':
            title = line.lstrip('>').strip()
            seq = ''
            while True:
                pos = f.tell()
                line = f.readline()
                if line == '':
                    yield title, seq
                    break
                elif line[0] == '>':
                    f.seek(pos)
                    yield title, seq
                    break
                else:
                    seq += line.strip()
        if line == '':
            break

def make_fasta_per_protein ():
    '''
    print(f'Canonical SARS-CoV2 protein names from GenBank:')
    for prot in canonical_protein_names:
        print(f'\t{prot}')
    '''
    f = open(args.input)
    records = {}
    for title, seq in fasta_parser(f):
        #if 'partial' in title:
        #    continue
        toks = title.split('|')
        gi = toks[0].strip()
        prot = toks[1].split('[')[0].strip().lower()
        if prot.startswith('chain '):
            continue # ignores PDB entries
        if prot in ['nsp11', 'structural protein']: # partial, ambiguous, etc. to ignore
            continue
        #if prot.endswith(' protein') and prot != 'structural protein':
        #    prot = prot[:-8]
        if prot.endswith(' proteiin'):
            prot = prot[:-9] + ' protein'
        #elif prot.endswith(' polyprotein'):
        #    prot = prot[:-12]
        elif prot.startswith('nonstructural protein '):
            prot = prot[22:]
        if prot in ['spike protein', 'spike glycoprotein', 's', 's protein']:
            prot = 'surface glycoprotein'
        elif prot in ['e', 'envelope', 'e protein']:
            prot = 'envelope protein'
        elif prot in ['m', 'membrane', 'membrane protein', 'matrix protein', 'matrix', 'm protein']:
            prot = 'membrane glycoprotein'
        elif prot in ['nucleocapsid protein', 'n', 'n protein']:
            prot = 'nucleocapsid phosphoprotein'
        elif prot in ['ns3', 'orf3', 'orf3 protein', 'orf3a']:
            prot = 'orf3a protein'
        elif prot in ['ns6']:
            prot = 'orf6 protein'
        elif prot in ['ns7a', 'orf7 protein']:
            prot = 'orf7a protein'
        elif prot in ['ns7b', 'orf7b protein']:
            prot = 'orf7b'
        elif prot in ['ns8']:
            prot = 'orf8 protein'
        #if prot in msa_grouping:
        #    msa_prot = msa_grouping[prot]
        #else:
        #    msa_prot = prot
        msa_prot = prot
        if msa_prot not in records:
            records[msa_prot] = {}
        records[msa_prot][f'{gi}:{prot}'] = seq
    prots = list(records.keys())
    prots.sort()
    noncanonical = []
    msa_targets = []
    num_seqs = {}
    for prot in prots:
        if prot not in canonical_protein_names:
            noncanonical.append(prot)
        giprots = list(records[prot].keys())
        if len(giprots) > 1:
            msa_targets.append(prot)
        else:
            print(f'# protein to be ignored due to single entry: {prot}')
    if len(noncanonical) > 0:
        print(f'{len(noncanonical)} non-canonical protein names found:')
        for prot in noncanonical:
            print(f'\t{prot}')
    if len(msa_targets) > 0:
        print(f'Creating FASTA files for multiple sequence alignment for proteins with multiple records:')
        for msa_target in msa_targets:
            giprots = list(records[msa_target].keys())
            giprots.sort()
            prot_nospace = msa_target.replace(' ', '_')
            fn = f'{prot_nospace}.fasta'
            print(f'\t{fn}')
            wf = open(fn, 'w')
            for giprot in giprots:
                wf.write(f'>{giprot}\n')
                wf.write(f'{records[msa_target][giprot]}\n')
            wf.close()
            num_seqs[msa_target] = len(giprots)
            #for gi in records[prot]:
            #    print(f'    {records[prot][gi]}')
    return msa_targets, num_seqs

async def get_clustered_proteins (request):
    f = open('clustered.txt')
    msa_targets = []
    for line in f:
        toks = line.strip().split('\t')
        msa_targets.append([toks[0].strip().replace(' ', '_'), int(toks[1])])
    f.close()
    return web.json_response(msa_targets)

def align ():
    msa_targets, num_seqs = make_fasta_per_protein()
    wf_l = open('clustered.txt', 'w')
    for msa_target in msa_targets:
        prot = msa_target.replace(' ', '_')
        if prot in args.skip_target:
            print(f'{prot} alignment is skipped by option.')
            continue
        cmd = ['.' + os.sep + 'clustalo', '-i', f'{prot}.fasta', '-o', f'{prot}.clustal', '--outfmt=clu', '--force']
        print(' '.join(cmd))
        pl = platform.platform()
        if pl.startswith('Windows'):
            subprocess.run(cmd, shell=True)
        elif pl.startswith('Linux'):
            subprocess.run(cmd)
        elif pl.startswith('Darwin'):
            subprocess.run(cmd)
        else:
            subprocess.run(cmd)
        wf_l.write(f'{msa_target}\t{num_seqs[msa_target]}\n')
    wf_l.close()

def on_shutdown (app):
    exit()

def view ():
    print(f'Starting web viewer for MSAs...')
    loop = asyncio.get_event_loop()
    app = web.Application()
    app.on_shutdown.append(on_shutdown)
    source_dir = os.path.dirname(os.path.realpath(__file__))
    app.router.add_route('GET', '/cll', get_clustered_proteins)
    app.router.add_static('/', os.path.join(source_dir))
    pl = platform.platform()
    if pl.startswith('Windows'):
        def_host = 'localhost'
    elif pl.startswith('Linux'):
        if 'Microsoft' in pl:
            def_host = 'localhost'
        else:
            def_host = '0.0.0.0'
    elif pl.startswith('Darwin'):
        def_host = '0.0.0.0'
    else:
        def_host = 'localhost'
    url = f'http://{def_host}:8080/index.html'
    loop.create_task(open_url(url))
    web.run_app(app)

async def open_url (url):
    webbrowser.open(url)

server_started = False
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input', default=None, help='SARS-CoV2 protein FASTA file, which can be downloaded at https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein&VirusLineage_ss=SARS-CoV-2,%%20taxid:2697049 (short URL: shorturl.at/tzIU3). Click "Download" button and choose "Sequence data (FASTA Format) Protein" at Step 1, "Download All Records" at Step 2, and "Use default: Accession GenBank Title" at Step 3. Then, click "Download" button. The downloaded file"s path should be given as an argument.')
parser.add_argument('--skip', dest='skip_target', nargs='*', default=[], help='If you want to skip alignment for certain proteins, write the proteins" names with this option. Protein names are all lower cases and space should be replaced with an underscore. For example, if you want to skip aligning ORF1ab polyprotein which takes the longest time, give --skip orf1ab_polyprotein. To see all canonical protein names, run python covid19.py --show-canonical')
parser.add_argument('--no-align', dest='no_align', action='store_true', default=False, help='Skips alignment')
parser.add_argument('--show-canonical', dest='show_canonical', action='store_true', default=False, help='Shows the canonical protein names which can be used in --skip option and quits.')
if sys.argv[0].startswith('python'):
    sys.argv = sys.argv[1:]
if sys.argv[0].endswith(os.path.basename(__file__)):
    sys.argv = sys.argv[1:]
args = parser.parse_args(sys.argv)

print(f'Canonical SARS-CoV2 protein names from GenBank:')
for prot in canonical_protein_names:
    print(f'\t{prot}')
if args.show_canonical:
    exit()
if args.no_align == False:
    if args.input is None:
        print(f'SARS-CoV2 protein FASTA file should be given with -i option to produce alignment, which can be downloaded at https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein&VirusLineage_ss=SARS-CoV-2,%%20taxid:2697049 (short URL: shorturl.at/tzIU3). Click "Download" button and choose "Sequence data (FASTA Format) Protein" at Step 1, "Download All Records" at Step 2, and "Use default: Accession GenBank Title" at Step 3. Then, click "Download" button. The downloaded file"s path should be given as an argument.')
        exit()
    else:
        align()
view()
