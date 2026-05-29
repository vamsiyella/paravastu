import sys
sys.path.insert(0, 'src')
from pipeline import run_batch

pdb_map = {
    17561: '2LBH',
    15409: '2KIB',
    16318: '2JWU',
    15380: '1LY2',
    15865: '1M8M',
    6838:  '1YMZ',
    17557: '2KSF',
    16299: '2JSV'
}

combined_stats = run_batch(
    bmrb_ids=list(pdb_map.keys()),
    pdb_map=pdb_map,
    run_dssp=True,
)