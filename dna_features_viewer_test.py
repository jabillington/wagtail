# %%

from dna_features_viewer import GraphicFeature, GraphicRecord
features=[
    GraphicFeature(start=0, end=20, strand=+1, color="#ffd700",
                   label="Small feature"),
    GraphicFeature(start=20, end=500, strand=+1, color="#ffcccc",
                   label="Gene 1 with a very long name"),
    GraphicFeature(start=400, end=700, strand=-1, color="#cffccc",
                   label="Gene 2"),
    GraphicFeature(start=600, end=900, strand=+1, color="#ccccff",
                   label="Gene 3")
]
record = GraphicRecord(sequence_length=1000, features=features)
record.plot(figure_width=5)

# %%




# %%


# %%
from geneblocks import load_record, DiffBlocks, CommonBlocks
import matplotlib.pyplot as plt

seq_1 = load_record("d378_attB-entry.gb")
seq_2 = load_record("d378_attB-mneon-puro-n.gb")

# FIND COMMON BLOCKS AND DIFFS
common_blocks = CommonBlocks.from_sequences({'seq_1': seq_1, 'seq_2': seq_2})
diff_blocks = DiffBlocks.from_sequences(seq_1, seq_2).merged()




# PLOT EVERYTHING
fig, axes = plt.subplots(3, 1, figsize=(15, 8))
common_blocks.plot_common_blocks(axes=axes[:-1])
diff_blocks.plot(ax=axes[-1], separate_axes=False)
axes[-1].set_xlabel("Changes in seq2 vs. seq1")
fig.savefig("complex_sequences.png", bbox_inches='tight')


common_blocks
# %%

# %%
