import argparse
import pandas as pd
from mygene import MyGeneInfo



def parse_args():
    parser = argparse.ArgumentParser(description="Perform analysis on the dataset.")
    parser.add_argument('--features-file-k562', type=str, required=True, help='Path to the input dataset (K562 features)')
    parser.add_argument('--features-file-hepg2', type=str, required=True, help='Path to the input dataset (HepG2 features)')
    parser.add_argument('--features-file-sknsh', type=str, required=True, help='Path to the input dataset (SK-N-SH features)')
    parser.add_argument('--rna-seq-file-k562', type=str, required=True, help='Path to the RNA-Seq data file')
    parser.add_argument('--rna-seq-file-hepg2', type=str, required=True, help='Path to the RNA-Seq data file')
    return parser.parse_args()


def main():
    args = parse_args()

    # Load the datasets
    features_df_k562 = pd.read_csv(args.features_file_k562, sep=r'\s+')
    features_df_hepg2 = pd.read_csv(args.features_file_hepg2, sep=r'\s+')
    features_df_sknsh = pd.read_csv(args.features_file_sknsh, sep=r'\s+')
    rna_seq_df_k562 = pd.read_csv(args.rna_seq_file_k562, sep='\t')
    rna_seq_df_hepg2 = pd.read_csv(args.rna_seq_file_hepg2, sep='\t')

    genes_k562 = features_df_k562['Feature'].tolist()
    mg = MyGeneInfo()
    gene_info = mg.querymany(genes_k562, scopes='symbol', fields='ensembl.gene', species='human')
    k562_feature_name_to_ensembl_id = {}
    for info in gene_info:
        try:
            gene_id = info['ensembl']['gene']
            k562_feature_name_to_ensembl_id[info['query']] = gene_id
        except TypeError:
            for entry in info['ensembl']:
                k562_feature_name_to_ensembl_id[info['query']] = entry['gene']
                break
        except KeyError:
            pass

    num_ensembl_ids_k562 = len([id for id in k562_feature_name_to_ensembl_id.values() if id is not None])
    print(num_ensembl_ids_k562, "Ensembl IDs retrieved for K562.")
    print(f"Could not retrieve Ensembl IDs for {len(genes_k562) - num_ensembl_ids_k562} features.")

    genes_hepg2 = features_df_hepg2['Feature'].tolist()
    gene_info = mg.querymany(genes_hepg2, scopes='symbol', fields='ensembl.gene', species='human')
    hepg2_feature_name_to_ensembl_id = {}
    for info in gene_info:
        try:
            gene_id = info['ensembl']['gene']
            hepg2_feature_name_to_ensembl_id[info['query']] = gene_id
        except TypeError:
            for entry in info['ensembl']:
                hepg2_feature_name_to_ensembl_id[info['query']] = entry['gene']
                break
        except KeyError:
            pass

    num_ensembl_ids_hepg2 = len([id for id in hepg2_feature_name_to_ensembl_id.values() if id is not None])
    print(num_ensembl_ids_hepg2, "Ensembl IDs retrieved for HepG2.")
    print(f"Could not retrieve Ensembl IDs for {len(genes_hepg2) - num_ensembl_ids_hepg2} features.")

    genes_sknsh = features_df_sknsh['Feature'].tolist()
    gene_info = mg.querymany(genes_sknsh, scopes='symbol', fields='ensembl.gene', species='human')
    sknsh_feature_name_to_ensembl_id = {}
    for info in gene_info:
        try:
            gene_id = info['ensembl']['gene']
            sknsh_feature_name_to_ensembl_id[info['query']] = gene_id
        except TypeError:
            for entry in info['ensembl']:
                sknsh_feature_name_to_ensembl_id[info['query']] = entry['gene']
                break
        except KeyError:
            pass

    num_ensembl_ids_sknsh = len([id for id in sknsh_feature_name_to_ensembl_id.values() if id is not None])
    print(num_ensembl_ids_sknsh, "Ensembl IDs retrieved for SK-N-SH.")
    print(f"Could not retrieve Ensembl IDs for {len(genes_sknsh) - num_ensembl_ids_sknsh} features.")

    gene_ids_to_tpm_k562 = {}
    for _, row in rna_seq_df_k562.iterrows():
        gene_id = row['gene_id']
        gene_id = gene_id.split('.')[0]
        tpm = row['TPM']
        gene_ids_to_tpm_k562[gene_id] = tpm

    gene_ids_to_tpm_hepg2 = {}
    for _, row in rna_seq_df_hepg2.iterrows():
        gene_id = row['gene_id']
        gene_id = gene_id.split('.')[0]
        tpm = row['TPM']
        gene_ids_to_tpm_hepg2[gene_id] = tpm

    k562_feature_to_k562_tpm = {}
    k562_feature_to_hepg2_tpm = {}

    for feature, ensembl_id in k562_feature_name_to_ensembl_id.items():
        if ensembl_id in gene_ids_to_tpm_k562:
            k562_feature_to_k562_tpm[feature] = gene_ids_to_tpm_k562[ensembl_id]
        else:
            k562_feature_to_k562_tpm[feature] = 0

        if ensembl_id in gene_ids_to_tpm_hepg2:
            k562_feature_to_hepg2_tpm[feature] = gene_ids_to_tpm_hepg2[ensembl_id]
        else:
            k562_feature_to_hepg2_tpm[feature] = 0

    f_out = open("k562_features_to_TPMs", "w")
    f_out.write("Feature\tTPM_K562\tTPM_HepG2\n")
    for feature, tpm in k562_feature_to_k562_tpm.items():
        tpm_hepg2 = k562_feature_to_hepg2_tpm[feature]
        tpm_k562 = k562_feature_to_k562_tpm[feature]
        f_out.write(f"{feature}\t{tpm_k562}\t{tpm_hepg2}\n")
    f_out.close()


    hepg2_feature_to_k562_tpm = {}
    hepg2_feature_to_hepg2_tpm = {}

    for feature, ensembl_id in hepg2_feature_name_to_ensembl_id.items():
        if ensembl_id in gene_ids_to_tpm_k562:
            hepg2_feature_to_k562_tpm[feature] = gene_ids_to_tpm_k562[ensembl_id]
        else:
            hepg2_feature_to_k562_tpm[feature] = 0

        if ensembl_id in gene_ids_to_tpm_hepg2:
            hepg2_feature_to_hepg2_tpm[feature] = gene_ids_to_tpm_hepg2[ensembl_id]
        else:
            hepg2_feature_to_hepg2_tpm[feature] = 0

    f_out = open("hepg2_features_to_TPMs", "w")
    f_out.write("Feature\tTPM_K562\tTPM_HepG2\n")
    for feature, tpm in hepg2_feature_to_k562_tpm.items():
        tpm_hepg2 = hepg2_feature_to_hepg2_tpm[feature]
        tpm_k562 = hepg2_feature_to_k562_tpm[feature]
        f_out.write(f"{feature}\t{tpm_k562}\t{tpm_hepg2}\n")
    f_out.close()

    sknsh_feature_to_k562_tpm = {}
    sknsh_feature_to_hepg2_tpm = {}

    for feature, ensembl_id in sknsh_feature_name_to_ensembl_id.items():
        if ensembl_id in gene_ids_to_tpm_k562:
            sknsh_feature_to_k562_tpm[feature] = gene_ids_to_tpm_k562[ensembl_id]
        else:
            sknsh_feature_to_k562_tpm[feature] = 0

        if ensembl_id in gene_ids_to_tpm_hepg2:
            sknsh_feature_to_hepg2_tpm[feature] = gene_ids_to_tpm_hepg2[ensembl_id]
        else:
            sknsh_feature_to_hepg2_tpm[feature] = 0

    f_out = open("sknsh_features_to_TPMs", "w")
    f_out.write("Feature\tTPM_K562\tTPM_HepG2\n")
    for feature, tpm in sknsh_feature_to_k562_tpm.items():
        tpm_hepg2 = sknsh_feature_to_hepg2_tpm[feature]
        tpm_k562 = sknsh_feature_to_k562_tpm[feature]
        f_out.write(f"{feature}\t{tpm_k562}\t{tpm_hepg2}\n")
    f_out.close()



    # now get features that are common in K562 and HepG2, and write out their TPMs in both cell lines
    common_features = set(features_df_k562['Feature']).intersection(set(features_df_hepg2['Feature']))
    f_out = open("features_in_k562_and_hepg2_TPMs", "w")
    f_out.write("Feature\tTPM_K562\tTPM_HepG2\n")
    for feature in common_features:
        ensembl_id = k562_feature_name_to_ensembl_id.get(feature)
        if ensembl_id is not None:
            tpm_k562 = gene_ids_to_tpm_k562.get(ensembl_id, 0)
            tpm_hepg2 = gene_ids_to_tpm_hepg2.get(ensembl_id, 0)
            f_out.write(f"{feature}\t{tpm_k562}\t{tpm_hepg2}\n")
    f_out.close()


    # now get features that are common in K562 and SKNSH, and write out their TPMs in both cell lines
    common_features = set(features_df_k562['Feature']).intersection(set(features_df_sknsh['Feature']))
    f_out = open("features_in_k562_and_sknsh_TPMs", "w")
    f_out.write("Feature\tTPM_K562\tTPM_HepG2\n")
    for feature in common_features:
        ensembl_id = k562_feature_name_to_ensembl_id.get(feature)
        if ensembl_id is not None:
            tpm_k562 = gene_ids_to_tpm_k562.get(ensembl_id, 0)
            tpm_hepg2 = gene_ids_to_tpm_hepg2.get(ensembl_id, 0)
            f_out.write(f"{feature}\t{tpm_k562}\t{tpm_hepg2}\n")
    f_out.close()

    # now get features that are common in HepG2 and SKNSH, and write out their TPMs in both cell lines
    common_features = set(features_df_hepg2['Feature']).intersection(set(features_df_sknsh['Feature']))
    f_out = open("features_in_hepg2_and_sknsh_TPMs", "w")
    f_out.write("Feature\tTPM_K562\tTPM_HepG2\n")
    for feature in common_features:
        ensembl_id = hepg2_feature_name_to_ensembl_id.get(feature)
        if ensembl_id is not None:
            tpm_k562 = gene_ids_to_tpm_k562.get(ensembl_id, 0)
            tpm_hepg2 = gene_ids_to_tpm_hepg2.get(ensembl_id, 0)
            f_out.write(f"{feature}\t{tpm_k562}\t{tpm_hepg2}\n")
    f_out.close()

    # now get features that are common in K562, HepG2, and SKNSH, and write out their TPMs in all 3 cell lines
    common_features = set(features_df_k562['Feature']).intersection(set(features_df_hepg2['Feature'])).intersection(set(features_df_sknsh['Feature']))
    f_out = open("features_in_k562_and_hepg2_and_sknsh_TPMs", "w")
    f_out.write("Feature\tTPM_K562\tTPM_HepG2\n")
    for feature in common_features:
        ensembl_id = k562_feature_name_to_ensembl_id.get(feature)
        if ensembl_id is not None:
            tpm_k562 = gene_ids_to_tpm_k562.get(ensembl_id, 0)
            tpm_hepg2 = gene_ids_to_tpm_hepg2.get(ensembl_id, 0)
            f_out.write(f"{feature}\t{tpm_k562}\t{tpm_hepg2}\n")
    f_out.close()

    # now get features that are unique to K562, and write out their TPMs in all 3 cell lines
    unique_k562_features = set(features_df_k562['Feature']).difference(set(features_df_hepg2['Feature'])).difference(set(features_df_sknsh['Feature']))
    f_out = open("features_unique_to_k562_TPMs", "w")
    f_out.write("Feature\tTPM_K562\tTPM_HepG2\n")
    for feature in unique_k562_features:
        ensembl_id = k562_feature_name_to_ensembl_id.get(feature)
        if ensembl_id is not None:
            tpm_k562 = gene_ids_to_tpm_k562.get(ensembl_id, 0)
            tpm_hepg2 = gene_ids_to_tpm_hepg2.get(ensembl_id, 0)
            f_out.write(f"{feature}\t{tpm_k562}\t{tpm_hepg2}\n")
    f_out.close()

    # now get features that are unique to HepG2, and write out their TPMs in all 3 cell lines
    unique_hepg2_features = set(features_df_hepg2['Feature']).difference(set(features_df_k562['Feature'])).difference(set(features_df_sknsh['Feature']))
    f_out = open("features_unique_to_hepg2_TPMs", "w")
    f_out.write("Feature\tTPM_K562\tTPM_HepG2\n")
    for feature in unique_hepg2_features:
        ensembl_id = hepg2_feature_name_to_ensembl_id.get(feature)
        if ensembl_id is not None:
            tpm_k562 = gene_ids_to_tpm_k562.get(ensembl_id, 0)
            tpm_hepg2 = gene_ids_to_tpm_hepg2.get(ensembl_id, 0)
            f_out.write(f"{feature}\t{tpm_k562}\t{tpm_hepg2}\n")
    f_out.close()

    # now get features that are unique to SKNSH, and write out their TPMs in all 3 cell lines
    unique_sknsh_features = set(features_df_sknsh['Feature']).difference(set(features_df_k562['Feature'])).difference(set(features_df_hepg2['Feature']))
    f_out = open("features_unique_to_sknsh_TPMs", "w")
    f_out.write("Feature\tTPM_K562\tTPM_HepG2\n")
    for feature in unique_sknsh_features:
        ensembl_id = sknsh_feature_name_to_ensembl_id.get(feature)
        if ensembl_id is not None:
            tpm_k562 = gene_ids_to_tpm_k562.get(ensembl_id, 0)
            tpm_hepg2 = gene_ids_to_tpm_hepg2.get(ensembl_id, 0)
            f_out.write(f"{feature}\t{tpm_k562}\t{tpm_hepg2}\n")
    f_out.close()

if __name__ == "__main__":    
    main()