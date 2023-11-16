import os
import json
import networkx as nx
import pandas as pd

from utils import annotate_cofactors, annotate_chemical_svg, get_autonomous_html, json_format, source_target, pathways_info_maker, insert_paths_ids
from Viewer import Viewer

if __name__ == '__main__':

    # Arguments
    d_idx = 0
    p_idx = 27
    pair_name = 'D'+str(d_idx)+'P'+str(p_idx)
    input_JSON = '/home/hector/Doctorado/MonteCarlo/Updated_graph/Graph_'+ pair_name +'.json'
    output_folder = pair_name
    template_folder = 'templates'
    autonomous_html = pair_name + '.html'
    # Make out folder if needed
    if not os.path.isfile(output_folder):
        try:
            os.makedirs(output_folder, exist_ok=True)
        except IOError as e:
            raise e

    data = json.load(open(input_JSON))
    G = nx.node_link_graph(data)
    inicial = pd.read_csv('/home/hector/Initial_met.csv',usecols=['SMILES'])
    inicial = inicial['SMILES'].tolist()
    producible = pd.read_csv('/home/hector/Doctorado/metabolitos_producibles.csv',usecols=['Smile'])
    producible = list(dict.fromkeys(producible['Smile']))
    producible.remove(producible[-1])

    all_paths = []
    for i in nx.shortest_simple_paths(G, producible[p_idx], inicial[d_idx]):
        all_paths.append(i)
    for i in range(len(all_paths)):
        for j in range(1,len(all_paths[i])):
            if '_' not in all_paths[i][j-1]:
                all_paths[i][j-1] = list(G.predecessors(all_paths[i][j]))
            else:
                all_paths[i][j-1] = [all_paths[i][j-1]]
        all_paths[i][-1] = [all_paths[i][-1]]
        
    network = json_format(data)
    
    pathways_info = pathways_info_maker(all_paths)
    network = insert_paths_ids(network, pathways_info)
    network = source_target(network, producible[p_idx], inicial[d_idx])
    # Add annotations
    # network = annotate_cofactors(network, args.cofactor)  # Typical cofactors
    network = annotate_chemical_svg(network)  # SVGs depiction for chemical

    # Build the Viewer
    viewer = Viewer(out_folder=output_folder, template_folder=template_folder)
    viewer.copy_templates()

    # Write info extracted from rpSBMLs
    json_out_file = os.path.join(output_folder, 'network.json')
    with open(json_out_file, 'w') as ofh:
        ofh.write('network = ' + json.dumps(network, indent=4))
        ofh.write(os.linesep)
        ofh.write('pathways_info = ' + json.dumps(pathways_info, indent=4))
    
    # Write single HTML if requested
    if autonomous_html is not None:
        str_html = get_autonomous_html(output_folder)
        with open(autonomous_html, 'wb') as ofh:
            ofh.write(str_html)
