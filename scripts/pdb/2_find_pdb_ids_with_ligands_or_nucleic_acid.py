import requests

from pdb_common import get_all_ligand_hetcodes
from pdb_common import get_pdb_ids_file_path

API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
API_PAGE_SIZE = 1000

def get_pdb_ids_page(start: int, ligands: list[str]):
    """Queries RCSB PDB API for proteins in complex with a nucleobase-containing small molecule or nucleic acid."""

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "label": "text",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": 2.5
                    }
                },
                {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "rcsb_entry_info.selected_polymer_entity_types",
                                "operator": "exact_match",
                                "value": "Protein/NA"
                            }
                        },
                        {
                            "type": "group",
                            "logical_operator": "and",
                            "nodes": [
                                {
                                    "type": "terminal",
                                    "service": "text",
                                    "parameters": {
                                        "attribute": "rcsb_entry_info.selected_polymer_entity_types",
                                        "operator": "exact_match",
                                        "value": "Protein (only)"
                                    }
                                },
                                {
                                    "type": "group",
                                    "logical_operator": "or",
                                    "nodes": [
                                        {
                                            "type": "group",
                                            "logical_operator": "or",
                                            "label": "nested-attribute",
                                            "nodes": [
                                                {
                                                    "type": "terminal",
                                                    "service": "text",
                                                    "parameters": {
                                                        "attribute": "rcsb_nonpolymer_instance_annotation.comp_id",
                                                        "operator": "exact_match",
                                                        "value": ligand
                                                    }
                                                }
                                            ]
                                        } for ligand in ligands
                                    ]
                                }
                            ]
                        }
                    ]
                },
            ]
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": start,
                "rows": API_PAGE_SIZE
            },
            "results_content_type": [
                "experimental"
            ]
        }
    }

    response = requests.post(API_URL, json=query)
    response.raise_for_status()

    results = response.json()
    pdb_ids = [hit["identifier"].lower() for hit in results.get("result_set", [])]
    total_count = results['total_count']
    next_page_start = start + API_PAGE_SIZE if start + API_PAGE_SIZE < total_count else None
    return pdb_ids, next_page_start


if __name__ == '__main__':
    ligands = get_all_ligand_hetcodes()
    all_pdb_ids = []
    start = 0
    while start is not None:
        pdb_ids_page, start = get_pdb_ids_page(start, ligands)
        all_pdb_ids.extend(pdb_ids_page)

    get_pdb_ids_file_path().write_text('\n'.join(all_pdb_ids))
