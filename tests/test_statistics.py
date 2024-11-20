###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import shutil
import tempfile
import unittest

import pandas as pd
from sciutil import *

u = SciUtil()
from scipy.stats import mannwhitneyu
from tqdm import tqdm

from levseq.variantcaller import *
from levseq.simulation import *
from levseq.utils import *
from sklearn.preprocessing import OneHotEncoder


class TestClass(unittest.TestCase):

    @classmethod
    def setup_class(self):
        local = True
        # Create a base object since it will be the same for all the tests
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))

        self.data_dir = os.path.join(THIS_DIR, 'test_data/')
        if local:
            self.tmp_dir = os.path.join(THIS_DIR, 'test_data/tmp/')
            if os.path.exists(self.tmp_dir):
                shutil.rmtree(self.tmp_dir)
            os.mkdir(self.tmp_dir)
        else:
            self.tmp_dir = tempfile.mkdtemp(prefix='test_data')

    @classmethod
    def teardown_class(self):
        shutil.rmtree(self.tmp_dir)


class TestStats(TestClass):

    def test_stats(self):
        # Read in once we have a bam file etc
        ref_str = 'T'
        well_df = get_reads_for_well('ref', 'test_data/brian/GCACTGGACAGACGTAGG.bam', ref_str,
                                     msa_path=f'msa.fa')
        well_df.to_csv('test_brian.csv')

    def test_stats_df(self):
        # Get them w.r.t to a mutation
        processed_plate_df = pd.DataFrame()
        parent = '#PARENT#'
        value_column = 'pdt'
        normalise = True

        parent_values = processed_plate_df[processed_plate_df['Mutations'] == parent][value_column].values
        parent_mean = np.mean(parent_values)
        parent_sd = np.std(parent_values)

        # if nomrliase normalize with standard normalisation
        if normalise:
            processed_plate_df[f'{value_column} standard norm'] = (processed_plate_df[
                                                                       value_column].values - parent_mean) / parent_sd
            value_column = f'{value_column} standard norm'

        parent_values = list(processed_plate_df[processed_plate_df['Mutations'] == parent][value_column].values)
        parent_mean = np.mean(parent_values)
        parent_sd = np.std(parent_values)
        sd_cutoff = 1.5  # The number of standard deviations we want above the parent values
        # Now for all the other mutations we want to look if they are significant, first we'll look at combinations and then individually
        grouped_by_mutations = processed_plate_df.groupby('Mutations')

        rows = []
        for mutation, grp in tqdm(grouped_by_mutations):
            # Get the values and then do a ranksum test
            if mutation != parent:
                vals = list(grp[value_column].values)
                U1, p = None, None
                # Now check if there are 3 otherwise we just do > X S.D over - won't be sig anyway.
                if len(grp) > 2:
                    # Do stats
                    U1, p = mannwhitneyu(parent_values, vals, method="exact")
                mean_vals = np.mean(vals)
                sig = mean_vals > ((sd_cutoff * parent_sd) + parent_mean)
                rows.append([mutation, len(grp), mean_vals, mean_vals - parent_mean, sig, U1, p])
        stats_df = pd.DataFrame(rows, columns=['mutation', 'number of wells with mutation', 'mean',
                                               'amount greater than parent mean', f'greater than > {sd_cutoff} parent',
                                               'man whitney U stat', 'p-value'])
        stats_df # Save to use in the next one...

    def test_mutation_summary(self):
        from collections import defaultdict
        stats_df = pd.DataFrame()

        mutation_dict = defaultdict(list)
        for mutation in stats_df['mutation'].values:
            mutations = mutation.split('_')
            for m in mutations:
                mutation_dict[m].append(mutation)

        rows = []
        with pd.ExcelWriter('mutations.xlsx', engine='xlsxwriter') as writer:
            for mutation, mutations in mutation_dict.items():
                # Here we want to now get the values for each of these i.e. the stats values for each one and summarize it maybe for now we'll just make a excel file
                df1 = stats_df[stats_df['mutation'].isin(mutations)]
                mutation = mutation.replace('*', '.')
                df1.to_excel(writer, sheet_name=mutation)
                # Also just take the mean of the mean lol and the sum of the number of the wells
                rows.append(
                    [mutation, np.sum(df1['number of wells with mutation'].values), '|'.join(set(list(mutations))),
                     np.mean(df1['mean'].values),
                     np.median(df1['median'].values), np.mean(df1['amount greater than parent mean'].values),
                     np.max(df1['amount greater than parent mean'].values)])

        df = pd.DataFrame(rows, columns=['mutation', 'number of wells with mutation',
                                         'other-mutations', 'mean', 'median',
                                         'mean amount greater than parent', 'max amount greater than parent'])

    def test_pca(self):
        processed_plate_df = pd.DataFrame()
        aas = ''
        seqs = []
        one_hots_nc = []
        one_hots_aa = []
        # Initialize OneHotEncoder
        encoder = OneHotEncoder()
        encoder.fit(np.array(['A', 'T', 'G', 'C', '-', '*']).reshape(-1, 1))

        encoder_aa = OneHotEncoder()
        encoder_aa.fit(np.array(aas).reshape(-1, 1))

        for nc in processed_plate_df['nc_variant'].values:
            seq = translate(nc)
            one_hot_encoded = encoder.transform(np.array(list(nc)).reshape(-1, 1))
            one_hot_encoded_array = one_hot_encoded.toarray().flatten()
            one_hots_nc.append(one_hot_encoded_array)

            one_hot_encoded = encoder_aa.transform(np.array(list(seq)).reshape(-1, 1))
            one_hot_encoded_array = one_hot_encoded.toarray().flatten()
            one_hots_aa.append(one_hot_encoded_array)

        from sklearn.decomposition import PCA
        import matplotlib.pyplot as plt
        import seaborn as sns
        pca = PCA(n_components=20)
        X = np.array(one_hots_nc)
        pca = pca.fit(X)
        pcs = pca.transform(X)

        processed_plate_df['PC 1'] = pcs[:, 0]
        processed_plate_df['PC 2'] = pcs[:, 1]

        sns.scatterplot(processed_plate_df, x='PC 1', y='PC 2', hue='pdt standard norm')
        PC_values = np.arange(pca.n_components_) + 1
        plt.plot(PC_values, (pca.explained_variance_ratio_ * 100), 'o-', linewidth=2, color='blue')
        plt.title('Scree Plot')
        plt.xlabel('Principal Component')
        plt.ylabel('Variance Explained')
        plt.show()