# import needed libraries
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import pickle
import random
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.decomposition import TruncatedSVD
from tqdm import tqdm


def similarity_search(matrix, index_node, top_n):
    """
    Function takes as input a matrix, an integer representing a node id, and an integer representing the
    number of similar nodes to return. The function uses this information and calculates the cosine similarity
    between the index node and all other included nodes. The results are sorted and returned as a list of
    lists where each list contains a node identifier and the cosine similarity score the top set of similar as
    indicated by the input argument are returned.
    :param matrix: where each row represents a node and each column represent the node embedding
    :param index_node: an integer representing a node id
    :param top_n: an integer representing the number of similar nodes to return
    :return:
        similar_nodes: a list of lists where each list contains a node identifier and the cosine similarity
        score the top set of similar as indicated by the input argument are returned
    """
    # http://markhneedham.com/blog/2016/07/27/scitkit-learn-tfidf-and-cosine-similarity-for-computer-science-papers/

    # calculate similarity
    cosine_similarities = cosine_similarity(matrix[index_node:index_node + 1], matrix).flatten()
    rel_node_indices = [i for i in cosine_similarities.argsort()[::-1] if i != index_node]
    similar_nodes = [(node, cosine_similarities[node]) for node in rel_node_indices][0:top_n]

    return similar_nodes
  
  
  def main():
    
    #######################################################################################
    # STEP 1: PROCESS KG NODE EMBEDDINGS
    #######################################################################################
    file_dir = 'analyses/preeclampsia/'
    # read in kg embeddings
    kg = pd.read_pickle(file_dir+ 'KG_Data_Embeddings')
    kg_labels = pd.read_csv(file_dir + 'KG_Node_Metadata.csv', header=0)
    kg = kg.merge(kg_labels[['id', 'label']], left_on='id', right_on='id', how='left').reset_index(drop=True)

   # create subsets
    kg_no_genes = kg[kg['grp'] != 'Genes'].reset_index(drop=True)
    kg_no_genes['grp2'] = ['Entity'] * len(kg_no_genes)
    kg_genes = kg[kg['grp'] == 'Genes'].reset_index(drop=True)
    kg_genes['grp2'] = ['Genes'] * len(kg_genes)

    # ignorome genes
    ignorome_genes = pd.read_csv(file_dir + 'pe_ignorome.csv', sep=',', header=0)
    ig_genes = ignorome_genes[ignorome_genes['IGNOROME'].notna()]['IGNOROME'].astype(int).astype(str).to_frame()
    kw_genes = ignorome_genes[ignorome_genes['KNOW'].notna()]['KNOW'].astype(int).astype(str).to_frame()
    kw_genes = ignorome_genes[ignorome_genes['BOTH'].notna()]['BOTH'].astype(int).astype(str).to_frame()
    updt_groups = ['PE' if kg.iloc[x]['id'] in set(ig_genes['IGNOROME']) else kg.iloc[x]['grp'] for x in range(0, len(kg))]

    # map embeddings to kg node identifiers
    ig_gene_embeddings = ig_genes.merge(kg, left_on='IGNOROME', right_on='id', how='left').reset_index(drop=True)
    ig_gene_embeddings['grp2'] = ['Ignorome'] * len(ig_genes)
    kw_gene_embeddings = kw_genes.merge(kg, left_on='KNOW', right_on='id', how='left').reset_index(drop=True)
    kw_gene_embeddings['grp2'] = ['Knowome'] * len(kw_genes)

    # update labels for kg node identifiers
    kg_no_genes = pd.concat([kg_no_genes, ig_gene_embeddings, kw_gene_embeddings])
    kg_no_tsne = 'ignorome_all_tsne.npy'
    kg_all_tsne = pd.concat([kg_genes, ig_gene_embeddings, kw_gene_embeddings])
    kg_all_genes = 'ignorome_pe_all_genes_tsne.npy'
    pe_genes = pd.concat([ig_gene_embeddings, kw_gene_embeddings])

    
    #######################################################################################
    # STEP 2: VISUALIZE KG NODE EMBEDDINGS USING T-SNE
    #######################################################################################

    # prep data for t-SNE
    matrix = [list(x) for x in kg_no_genes['embeds']]
    matrix = [list(x) for x in kg_all_genes['embeds']]
    matrix = [list(x) for x in pe_genes['embeds']]
    full_matrix = [list(x) for x in kg['embeds']]

    # create t-sne
    X_reduced = PCA(n_components=25, random_state=1).fit_transform(matrix)
    X_embedded = TSNE(n_components=2, random_state=1, verbose=True, perplexity=50.0).fit_transform(X_reduced)

    # set up colors and legend labels
    names = {'Edge': 'Edge', 'Gene Ontology': 'Gene Ontology', 'Drugs': 'Drugs', 'Diseases': 'Diseases', 'Pathways': 'Pathways',
             'Genes': 'Genes', 'Phenotypes': 'Phenotypes', 'PE': 'PE Ignorome'}
    colors = {'Diseases': 'paleturquoise', 'Drugs': 'lavenderblush', 'Gene Ontology': 'palegreen', 'Genes': 'lightgrey',
              'Pathways': 'mistyrose', 'Phenotypes': 'lavender', 'PE Ignorome': 'cornsilk'}
    dis = mpatches.Patch(color='deepskyblue', label='Disease'); drg = mpatches.Patch(color='pink', label='Drug')
    go = mpatches.Patch(color='lightgreen', label='GO'); ge = mpatches.Patch(color='dimgray', label='Gene')
    pat = mpatches.Patch(color='crimson', label='Pathway'); phe = mpatches.Patch(color='purple', label='Phenotype')
    pe = mpatches.Patch(color='goldenrod', label='Ignorome')
    df = pd.DataFrame(dict(x=X_embedded[:, 0], y=X_embedded[:, 1], group=list(updt_groups)))
    groups = df.groupby('group')
    fig, ax = plt.subplots(figsize=(13, 10))
    for x, grp in groups:
        if x == 'Gene Ontology':
            ax.plot(grp.x, grp.y, marker='o', linestyle='', ms=4, label=names[x], color=colors[x], mec='forestgreen', alpha=0.6)
        if x == 'Genes':
            ax.plot(grp.x, grp.y, marker='o', linestyle='', ms=4, label=names[x], color=colors[x], mec='dimgray', alpha=0.6)
        if x == 'Drugs':
            ax.plot(grp.x, grp.y, marker='o', linestyle='', ms=4, label=names[x], color=colors[x], mec='hotpink', alpha=0.6)
        if x == 'Phenotypes':
            ax.plot(grp.x, grp.y, marker='o', linestyle='', ms=4, label=names[x], color=colors[x], mec='purple', alpha=0.6)
        if x == 'Pathways':
            ax.plot(grp.x, grp.y, marker='o', linestyle='', ms=4, label=names[x], color=colors[x], mec='crimson', alpha=0.6)
        if x == 'Diseases':
            ax.plot(grp.x, grp.y, marker='o', linestyle='', ms=4, label=names[x], color=colors[x], mec='deepskyblue', alpha=0.6)
        if x == 'PE':
            ax.plot(grp.x, grp.y, marker='*', linestyle='', ms=12, label=names['PE'], color=colors['PE Ignorome'], mec='darkgoldenrod', alpha=0.6)
    plt.legend(handles=[dis, drg, go, ge, pat, phe, pe], fontsize=12, frameon=False, loc="lower center", ncol=7)
    ax.tick_params(labelsize=13)
    m0 = -55; m1 = 55; plt.ylim(m0, m1); plt.xlim(m0, m1)
    plt.show(); plt.close()

    
    #######################################################################################
    # STEP 3: CALCULATE KG NODE EMBEDDING SIMILARITY
    #######################################################################################

    # create dicts to store information about each node and node type
    embedding_idx = kg['id'].to_dict(); embedding_id_grp = kg['grp'].to_dict()
    pe_gene_idx = {k: v for k, v in embedding_idx.items() if v in set(ig_genes['IGNOROME'])}
    kg_node_labels = {row['id']: row['label'] for idx, row in kg.iterrows() if row['grp'] != 'Edge'}

    # convert embeddings to compressed sparse matrix for cosine similarity calculation
    kg_matrix = csr_matrix([list(x) for x in kg['embeds']])

    # find similar entities
    similar_entities = {}
    for gene_idx, gene_id in tqdm(pe_gene_idx.items()):
        hit_dict = {}; gene_label = kg_node_labels[gene_id]
        matches = [x for x in similarity_search(kg_matrix, gene_idx, 200) if x[0] not in pe_gene_idx.keys()][0:100]
        if len(matches) < 100: break
        for x, score in matches:
            if embedding_id_grp[x] != 'Edge':
                match_id = embedding_idx[x]; node_grp = embedding_id_grp[x]; match_label = kg_node_labels[match_id]
                if node_grp in hit_dict.keys(): hit_dict[node_grp] += [[match_id, match_label, round(score, 3)]]
                else: hit_dict[node_grp] = [[match_id, match_label, round(score, 3)]]
        hit_dict['gene_label'] = gene_label
        similar_entities[gene_id] = hit_dict

    # create lists of the top annotated entities by node type
    annotation_counts = {'Drugs': {}, 'Diseases': {}, 'Genes': {}, 'Pathways': {}, 'Phenotypes': {}, 'Gene Ontology': {}}
    master_list = dict()
    for k, v in tqdm(similar_entities.items()):
        gene_label = v['gene_label']; val = {x: y for x, y in v.items() if x != 'gene_label'}
        master_list[k] = {'Drugs': [], 'Diseases': [], 'Genes': [], 'Pathways': [], 'Phenotypes': [], 'Gene Ontology': []}
        for x in val.keys():
            for i in val[x]:
                if i[0] in annotation_counts[x].keys():
                    annotation_counts[x][i[0]]['scores'].append(i[-1])
                    annotation_counts[x][i[0]]['counts'] += 1
                    annotation_counts[x][i[0]]['genes'].append('{} ({})'.format(k, gene_label))
                else:
                    annotation_counts[x][i[0]] = {
                        'scores': [i[-1]], 'counts': 1, 'genes': ['{} ({})'.format(k, gene_label)]}
                master_list[k][x].append(i)

    # print top annotated entities
    for k, v in annotation_counts.items():
        for x in ['Drugs', 'Diseases', 'Genes', 'Gene Ontology', 'Pathways', 'Phenotypes']:
            k = x; v = annotation_counts[k]; overall_scores = [i for j in [v[x]['scores'] for x in v.keys()] for i in j]
            hs = sorted(list(set(overall_scores)), reverse=True)
            y = [((x, kg_node_labels[x]), v[x]['scores']) for x in v.keys()
                 if len([g for g in v[x]['scores'] if g in hs]) > 0]
            with open(file_dir + 'pkt_validation/pe_pkt_most_similar_{}.txt'.format(x), 'w') as f:
                for z in y:
                    hits = ', '.join(v[z[0][0]]['genes'])
                    try: f.write(str(z[0][0]) + '*' + str(z[0][1]) + '*' + str(z[1][0]) + '*' + hits + '\n')
                    except AttributeError: pass
            # printing most frequent by count
            count_list = [v[x]['counts'] for x in v.keys()]
            entity_counts = sorted(count_list, reverse=True); min_cnt = entity_counts[-1]; max_cnt = entity_counts[0]
            cut_off = [max_cnt] if entity_counts.count(max_cnt) > 1 else entity_counts
            he = [((x, kg_node_labels[x]), v[x]['counts']) for x in v.keys() if v[x]['counts'] in entity_counts]
            with open(file_dir + 'pkt_validation/pe_pkt_most_frequent_{}.txt'.format(x), 'w') as f:
                for u in he:
                    try: f.write(str(u[0][0]) + '*' + str(u[0][1]) + '*' + str(u[1]) + '\n')
                    except AttributeError: pass
 

    #######################################################################################
    # STEP 4: Perfoming Modified Enrichment Analysis and Statistical Testing
    #######################################################################################
    
    write_loc = 'analyses/preeclampsia/'
    gold_standard_annotations = {
        'Diseases': {'DOID_9870', 'DOID_1700', 'DOID_681', 'DOID_332', 'DOID_1094', 'DOID_5158', 'DOID_1926', 'DOID_9452'},
        'Drugs': {'D005486', 'C445526', 'C010634', 'D004962', 'C059714', 'D006003', 'C098288', 'D018927', 'D008277',
                  'D000068900'},
        'Genes': {'5352', '27094', '10576', '4709', '64983', '3047', '5931', '148867', '204962', '58480'},
        'Gene Ontology': {'GO_0070125', 'GO_0006882', 'GO_0005747', 'GO_0031462', 'GO_0005833', 'GO_0005516', 'GO_0007062',
                          'GO_0042393', 'GO_0000398', 'GO_0051015'},
        'Pathways': {'R-HSA-6799198', 'R-HSA-5419276', 'R-HSA-611105', 'R-HSA-1566948', 'R-HSA-3906995', 'R-HSA-194840',
                     'R-HSA-391251', 'R-HSA-2500257', 'R-HSA-913531', 'R-HSA-212165'},
        'Phenotypes': {'HP_0008316', 'HP_0002725', 'HP_0008344', 'HP_0011904', 'HP_0001935'}
    }

    # read kg files
    kg = pd.read_pickle(file_dir + 'KG_Data_Embeddings')
    kg_labels = pd.read_csv(file_dir + 'KG_Node_Metadata.csv', header=0)
    kg = kg.merge(kg_labels[['id', 'label']], left_on='id', right_on='id', how='left').reset_index(drop=True)
    
    # set-up gene lists
    ig_genes = ignorome_genes[ignorome_genes['IGNOROME'].notna()]['IGNOROME'].astype(int).astype(str).to_frame()
    kw_genes = ignorome_genes[ignorome_genes['KNOW'].notna()]['KNOW'].astype(int).astype(str).to_frame()
    both_genes = ignorome_genes[ignorome_genes['BOTH'].notna()]['BOTH'].astype(int).astype(str).to_frame()
    pe_gene_idx = set(ig_genes['IGNOROME']) | set(kw_genes['KNOW']) | set(both_genes['BOTH'])
    create dicts to help with using matrix indices
    embedding_idx = kg['id'].to_dict(); embedding_id_grp = kg['grp'].to_dict()
    rev_embedding_idx = {v: k for k, v in embedding_idx.items()}
    kg_node_labels = {row['id']: row['label'] for idx, row in kg.iterrows() if row['grp'] != 'Edge'}
    lab_dicts = {'embedding_idx': embedding_idx, 'rev_embedding_idx': rev_embedding_idx, 'kg_node_labels' : kg_node_labels,
                 'embedding_id_grp': embedding_id_grp}
    embedding_idx = lab_dict['embedding_idx']; rev_embedding_idx = lab_dict['rev_embedding_idx']
    kg_node_labels = lab_dict['kg_node_labels']; embedding_id_grp = lab_dict['embedding_id_grp']

    # convert embeddings to compressed sparse matrix
    kg_matrix = csr_matrix([list(x) for x in kg['embeds']])

    # pull 1,000 samples (with replacement) from all genes other than PE
    sample_set = set(kg[kg['grp'] == 'Genes']['id']) - pe_gene_idx  # 22385
    genes = [random.choices(list(sample_set), k=445) for _ in range(0, 1000)]
    gene_samples = set([i for u in genes for i in u])
    updated = set([x for x in gene_samples if x not in gene_samples])
    sample_start = 0; sample_end = 101; genes = gene_samples[sample_start:sample_end]

    # find similar entities
    idx = list(range(0, 101))
    for samp in tqdm(range(sample_start, sample_end)):
        sampled_results = {}
        for gene_id in tqdm(genes[idx.pop(0)]):
            print(gene_id, len(genes))
            gene_idx = rev_embedding_idx[gene_id]; hit_dict = {}; gene_label = kg_node_labels[gene_id]
            matches = similarity_search(kg_matrix, gene_idx, 100); sampled_results[samp] = {}
            for x, score in matches:
                if embedding_id_grp[x] != 'Edge':
                    match_id = embedding_idx[x]; node_grp = embedding_id_grp[x]; match_label = kg_node_labels[match_id]
                    if node_grp in hit_dict.keys(): hit_dict[node_grp] += [[match_id, match_label, round(score, 3)]]
                    else: hit_dict[node_grp] = [[match_id, match_label, round(score, 3)]]
            hit_dict['gene_label'] = gene_label
            sampled_results[samp][gene_id] = hit_dict

    # reads in the annotation runs and organizes them by domain --> hits --> entity-scores
    from glob import glob
    f = glob(file_dir + 'sample_sets/*sample*')[:-1]
    sampled_results = {k: v for d in [pickle.load(open(x, 'rb')) for x in tqdm(f)] for k, v in d.items()}
    annotation_counts = {'Drugs': {}, 'Diseases': {}, 'Genes': {}, 'Pathways': {}, 'Phenotypes': {}, 'Gene Ontology': {}}
    
    for s in tqdm(list(sampled_results.keys())[0:1000]):
        annotation_counts['Drugs'][s] = []; annotation_counts['Diseases'][s] = []
        annotation_counts['Genes'][s] = []; annotation_counts['Gene Ontology'][s] = []
        annotation_counts['Pathways'][s] = []; annotation_counts['Phenotypes'][s] = []
        for k, v in sampled_results[s].items():
            for x in v.keys():
                if x != 'gene_label':
                    for i in v[x]:
                        annotation_counts[x][s] += ['{}:{}'.format(i[0], i[2])]

    # for any repeating annotations within the same run, this method condenses them by the max similarity score
    hits1 = {'Drugs': [], 'Diseases': [], 'Genes': [], 'Pathways': [], 'Phenotypes': [], 'Gene Ontology': []}
    for k, v in tqdm(annotation_counts.items()):
        hits1[k] = {}; ent_dict = {}
        for run, out in v.items():
            for i in out:
                term, score = i.split(':')
                if term in ent_dict.keys(): ent_dict[term] = max([ent_dict[term], float(score)])
                else: ent_dict[term] = float(score)
        hits1[k] = ['{}:{}'.format(k, v) for k, v in ent_dict.items()]

    # gets final hits and scores tests
    hits = {'Drugs': [], 'Diseases': [], 'Genes': [], 'Pathways': [], 'Phenotypes': [], 'Gene Ontology': []}
    for k, v in tqdm(hits1.items()):
        gs = gold_standard_annotations[k]; entity_list, entity_scores = [], []
        tup_list = [tuple(x.split(':')) for x in v]
        sorted_tuples = sorted(tup_list, key=lambda x: x[1], reverse=True)
        hits_to_search = [x[0] for x in sorted_tuples]; hit_scores = [float(x[1]) for x in sorted_tuples]
        hit = list(gs.intersection(set(hits_to_search)))
        if len(hit) > 0: hits[k] += ['{}: {}'.format(x, kg_node_labels[x].lower()) for x in hit]

    for k, v in hits.items():
        print('{}: {}, p-value={}'.format(k.upper(), len(v), len(v)/float(1000)))
