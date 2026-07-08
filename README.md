<h1>Adipose Tissue–Associated Lipid Serum Signatures Modulate Efficacy and Toxicity of CAR-T Cell Therapy in Lymphoid Malignancies</h1>

<p><strong>R</strong></p>

<p>
  This repository contains the R code for the lipidomics–clinical outcomes project:
  <strong>development of serum lipid signatures associated with efficacy and toxicity following BCMA-directed CAR T-cell therapy in lymphoid malignancies</strong>.
  The workflow spans data wrangling, MS preprocessing (lamivudine-normalized relative abundance), feature selection, correlation with adipose compartments, outcomes modeling (CRS, PFS, OS), and figure generation, with optional GO enrichment from co-culture scRNA-seq.
</p>

<p><strong>Note:</strong> Example/synthetic input data and full scripts will be aligned with the published manuscript. Until then, only structure and non-proprietary utilities are included.</p>

<hr />

<h2>Requirements</h2>
<ul>
  <li><strong>R ≥ 4.3.0</strong></li>
</ul>

<h3>Key R packages</h3>

<p>
  <strong>I/O &amp; wrangling</strong><br />
  <code>readxl</code>, <code>openxlsx</code>, <code>tidyr</code>, <code>dplyr</code>, <code>purrr</code>, <code>stringr</code>, <code>reshape2</code>, <code>data.table</code> (optional), <code>RJSONIO</code>
</p>

<p>
  <strong>Visualization</strong><br />
  <code>ggplot2</code>, <code>ggpubr</code>, <code>ggrepel</code>, <code>ggprism</code>, <code>hrbrthemes</code>, <code>gridExtra</code>, <code>ggpattern</code>, <code>plotly</code>, <code>ComplexHeatmap</code>, <code>ggvenn</code>, <code>ggVennDiagram</code>, <code>RColorBrewer</code>
</p>

<p>
  <strong>Metabolomics &amp; feature selection</strong><br />
  <code>MetaboAnalystR</code>, <code>pls</code>, <code>igraph</code>
</p>

<p>
  <strong>Survival &amp; competing risks</strong><br />
  <code>survival</code>, <code>survminer</code>, <code>ggsurvfit</code>, <code>tidycmprsk</code>
</p>

<p>
  <strong>Classification / ROC / ML</strong><br />
  <code>pROC</code>, <code>randomForest</code>, <code>caret</code>
</p>

<p>
  <strong>Meta-analysis &amp; utilities</strong><br />
  <code>metafor</code>, <code>datasets</code>
</p>

<p><em>If some packages are not required for your run, comment them out in <code>0_packages_setup.R</code>.</em></p>

<h2>Workflow Summary</h2>
<ol>
  <li><strong>QC &amp; normalization</strong>: Remove failed features/samples by internal QC; normalize to <strong>lamivudine</strong>; export relative abundance (a.u.).</li>
  <li><strong>Feature selection</strong>: FC filtering and <strong>PLS-DA</strong> (VIP threshold configurable).</li>
  <li><strong>Correlation filtering</strong>: Retain metabolites associated with AT compartments (threshold configurable in <code>config/parameters.yml</code>).</li>
  <li><strong>Outcomes modeling</strong>: CRS (logistic), PFS/OS (Cox) with prespecified covariates; generate model summaries and forest plots.</li>
  <li><strong>Visualization</strong>: Volcano plots, longitudinal boxplots (median + 95% CI; points = individuals; lines = median over time), PLS-DA.</li>
</ol>

<h2>Contact</h2>
<p>
  <strong>David M. Cordas dos Santos</strong><br />
  &lt;<a href="mailto:davidm_cordasdossantos@dfci.harvard.edu">davidm_cordasdossantos@dfci.harvard.edu</a>&gt;
</p>
