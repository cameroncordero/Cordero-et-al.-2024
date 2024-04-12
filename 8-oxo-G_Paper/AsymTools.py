
from asymtools.annot import *
from asymtools.plotting import *

# Load and annotate maf file
# An example lung cancer maf can be downloaded from https://gdc.cancer.gov/about-data/publications/luad_2014
m = preprocess_maf('/media/cam/Data9/8-oxodG_Final_Analysis/maf_files/KBr_Treated.maf')

# Plot asymmetries of mutation counts
# twin_bar_txplot(m)
# twin_bar_repplot(m)

# Plot with correction of genomic content in mutations/Mb
# Choose 'exome' or 'genome' as appropriate
twin_bar_txplot(m,normalization='genome')
twin_bar_repplot(m,normalization='genome')

# In addition to the matplotlib axes, the plotting functions can also return a dataframe with the plotted values
# The columns summarize for each mutation type the total counts (n1,n2), normalized rates if using (r1,r2), and ratio of complementary mutations (ratio)
ax,res = twin_bar_repplot(m,normalization='genome')
print(res[1])


plt.show()
