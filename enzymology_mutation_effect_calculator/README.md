## Enzymology Mutation Performance Calculator


### Short Description

This tool estimates the effect of individual amino acid mutations for a given protein mutagenesis
library and associated assay data.


### Tool Description

This tool estimates the first-order effect of individual amino acid mutations for a given protein
mutagenesis library and associated phenotypic data. Estimation is performed by fitting a regularized
linear model, namely
[_ElasticNet_](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLS.fit_regularized.html),
to the provided data. Note that only a single phenotype can be evaluated at once. Regularization
parameters are auto-tuned using randomized grouped 5-fold cross validation.

This analysis assumes that each unique mutation is observed in combinations with other mutations
such that an average effect can be determined.  Thus, the tool is probably most useful for
evaluating deep random mutagenic libraries, or combinatorial consolidation libraries. That said,
the tool may be applied to saturation libraries as well.

**As of June 10, 2021 mutation positions are reported using 1-based indexing.**

#### Approach

For this analysis, we model sequence performance using the linear model of the form:

    y(seq) = const + m_1 * mutation_1(seq) + m_2 * mutation_2(seq) +  ...

where `const` is a constant and `mutation_i(seq)` is 0 or 1 depending on whether `mutation_i`
e.g. `H5k`, is present in the given sequence `seq`, and the coefficients `m_i` are the individual
mutation effects.

Mutation interaction effects may also be modeled by changing the `degree` setting in the pipeline
parameters json file from its default value of `1` to `2`.  In this case, the following terms are
additionally included in the model:

    y(seq) += m_1_2 * mutation_1_2(seq) + m_1_3 * mutation_1_3(seq) + ...

where `mutation_i_j(seq)` is 0 or 1 depending on whether `mutation_i` and `mutation_j` e.g. `H5K`
and `G10V`, are both present in the given sequence `seq`, and the coefficients `m_i_j` are the
pairwise mutation effects. Note that higher order terms in the model are dropped if they are
redundant with lower order terms.

The mutation effect sizes are determined according to the following procedure:

1. Align mutant sequences to the reference sequence.
2. Identify point mutations in each of the mutant sequences, relative to the reference.
3. Create feature vectors for each sequence by performing 1-hot encoding of sequence mutations.
4. Determine optimal L1 and L2 regularization parameters for the linear model by shuffled grouped
   5-fold cross validation, where grouping is done by strain id.
5. Refit the model over the full data-set, using the tuned regularization parameters.
6. Report on model fit and individual mutation effects.


### Inputs

+  **Sequence Performance File (csv)**: A CSV file containing assay data. The first column should
   contain sequence identifiers, such as Strain ID, DnaComponent ID, or Sample ID.  At least one
   other column should contain the performance data for the analysis. It is generally assumed,
   though not strictly required, that performance data consists of fold-improvement-over-parent
   (FIOP) values. _Note that if you are running the tool using the `auto` sequence input mode, you
   must use the Strain ID as the sequence identifier for now._

+  **Assay Data Column Name**: The name of the column containing the performance data.

+  **Sequence Input Mode**: Select either `auto` or `manual`. In `manual` mode, you must provide
   sequence data directly using the following input:

   -  **Sequences File (fasta)** : A fasta file containing sequences of all proteins considered in
      the analysis. Sequences may either be expressed in terms of DNA or amino acids. Sequence IDs
      should match those in the Sequence Performance File.  _Note that this tool assumes that the
      first sequence in the file is the reference sequence, against which all mutant sequences will
      be evaluated._

   In `auto` mode, the performance file sequence identifiers are used to query mutated sequences
   directly from LIMS via the Genetic Changes Table. In addition, you can specify the reference
   sequence against which the mutant sequence should be evaluated using the following inputs:

   -  **Reference Sequence Id**: The Strain ID or DnaComponent ID associated with the reference
      sequence, against which mutant sequence will be compared.  If the ID does not correspond
      directly to the sequence of interest, you can you the subsequent options to sub-select from
      this sequence. The input also accepts two magic _string_ values, available from the input's
      drop down menu. These are `use_logical_parent_strain` and `use_reference_strain`.  The former
      attempts to find a common logical parent strain for all mutant strains. The later attempts to
      find a common earliest ancestral strain for all mutant strains.

   -  **Feature Term**: Optionally, filter sequence data to only include the sub-sequence with
      annotation of this type.

   -  **Annotation Name**: Optionally, filter sequence data to only include the sub-sequence with
      annotation matching this name.

   All queried sequences are saved to the input directory using a file name derived from the
   sequence performance file input.

+  **Sequence Translation Mode**: Select either `standard` or `custom`. In `standard` mode, the
   "standard" NCBI codon translation table is used to translate DNA inputs to their corresponding
   protein
   sequences. See [The Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for
   additional information.  In `cutsom` mode, you must provide a custom codon translation table
   using the following input:.

   -  **Codon Table (string)**: Provide a custom codon translation table using the following format:

              AAs    = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
              Starts = ---M------**--*----M---------------M----------------------------
              Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
              Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
              Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

      Note that the characters used in the `AAs` row need not adhere to the standard protein
      alphabet.

+ **Modeling Pipeline Parameters (json)**: A json formatted file containing additional parameters
  used by the modeling pipeline.  A pre-populated __pipeline_parameters.json__ file is provided in
  the __user_inputs__ directory and can typically be used as-is.  In certain rare-cases it may be
  useful to modify these parameters to accommodate especially noisy data. Please reach out to the
  Data Science team for additional information on these advanced settings.

+ **Library Type**: Select either `Single Site Library (SSL)` or `Combinatorial Library`.  When
  `SSL` is selected, results specific to Combinatorial Libraries, e.g. the Mutation Co-Occurance
  plot, are disabled.

+ **Coefficient Cutoff Probability**: Remove fitted mutation coefficients from the mutation
  coefficients figure if the associated t-statistic exceeds this probability threshold.  Note that
  all coefficients values are saved to the __model_coefficients__ output file, regardless of this
  setting.


### Outputs

This tool creates the following primary outputs:

+  **Model Summary**: A table of summary metrics, describing the model's overall fit. This
   includes quantities such as R^2, tuned hyper-parameter values, etc.

+  **Model Fit Plot**: A scatter plot of the actual sequence performance vs the predicted sequence
   performance.

+  **Model Coefficients**: A table of estimated mutation effects, significance levels,
   and confidence intervals. Note that only "significant" effects (> 4 standard deviations) are
   displayed.  The complete list of values is available in the _output\_files_ directory.

+  **Model Coefficients Plot**: A bar chart of signification mutation effects, arranged from most
   beneficial to least beneficial.

In addition, the tool creates a few diagnostic figures:

+  **Library Coverage**: A heatmap mutation frequency within the library, by position
   and type.

+  **Mutation Co-Occurrence**: A heatmap of pairwise mutation co-occurrence within the library.

+  **Sequence Distance**: A heatmap of pairwise sequence distances within the library, measured
   using the hamming metric.

Finally, a few supplemental artifacts are written to _output\_files_ directory:

+  **Model**: A pickled version of the fitted model.

+  **User Inputs**: A json file of user input values.

+  **Alignment (ali)**: The multiple-sequence alignment of all translated input sequences, using
   the clustal format.

+  **Alignment (csv)**: A table summarizing individual  mutations identified by the
   multiple-sequence alignment.

+  **Codon Table (txt)**: A codon translation table, if a custom codon table was used for
   translation.


If you are satisfied with the analysis, you can click the **Save Results To LIMS** button and the
contents of both the _input_ and _output_ file directories are copied to a new LIMS Dataset.


### Important Notes and Limitations

It is well known that protein stability and function are enabled by complex interactions between
proximal and distant amino acids.  While a linear model can identify potentially useful "average"
effects of specific mutations, this is by no means a guarantee that the combinations of such
mutations will be additive.  Nor should it be assumed that non-linear effects present in the data
_but not captured by the linear model_ can't provide greater overall benefits.  Thus, results
produced by this tool should be used only in conjunction with other optimization strategies that do
exploit (implicitly or explicitly) the underlying non-linearity in the data.


### Developer Notes

A few test data sets are located in the test directory.  Specifically, these are:

+  **.code/tests/fixtures/random\_proteins\_simulated**: a simulated library containing one-thousand
   random 5 change mutants of a 100 amino acid protein sequence, randomly seeded with 25(?)
   beneficial mutations. Run-time is ~5 minutes.
+  **.code/tests/fixtures/random\_proteins\_simulated\_small**: a simulated library containing
   one-thousand random 3 change mutants of a 10 amino acid protein sequence, randomly seeded with 5
   beneficial mutations. Run-time is ~30 seconds.
+  **.code/tests/fixtures/omicron\_avie**: a real-life combinatorial library created for the
   Omicron-Avie program.  It contains around 800 mutants, obtained by randomly combining ~25
   beneficial mutations previous selected from a saturation library into the wildtype sequence.
   Run-time is ~5 minutes.
