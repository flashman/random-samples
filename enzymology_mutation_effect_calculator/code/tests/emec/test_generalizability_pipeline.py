import unittest

from emec import pipeline
from emec import generalizability_pipeline
from emec.utils.simulation_utils import RandomMutagenesisSimulation


class TestGeneralizabilityPipeline(unittest.TestCase):
    def test_run_pipline(self):
        # For this test, ensure that beneficial mutations are always well separated from neutral
        # mutations.

        # Simulate data.
        mutation_params = {"n_mutations": 3}
        fiop_params = {"n_beneficial": 5, "benefit_mu": 2.0, "benefit_sigma": 0.05}

        sim = RandomMutagenesisSimulation(10, 999, mutation_params, fiop_params)
        sim.simulate_sequences()
        sim.simulate_fiops_v1()

        # Run the evaluation pipeline
        results, fitted_params, elastic_net = pipeline.run_pipeline(
            sim.alignment, sim.attributes_df["fiop"], pipeline.PIPELINE_PARAMS
        )

        params = generalizability_pipeline.GeneralizabilityParams(n_splits=2, n_steps=4)
        (
            model_fit_df,
            model_coef_df,
            feature_counts_df,
        ) = generalizability_pipeline.run_pipeline(
            results, params, pipeline.PIPELINE_PARAMS
        )

        self.assertEqual(len(model_fit_df), 8)
        self.assertEqual(len(model_coef_df), 8)
        self.assertEqual(len(feature_counts_df), 8)

    def test_plots(self):
        # Simulate data.
        mutation_params = {"n_mutations": 3}
        fiop_params = {"n_beneficial": 5, "benefit_mu": 2.0, "benefit_sigma": 0.05}

        sim = RandomMutagenesisSimulation(8, 200, mutation_params, fiop_params)
        sim.simulate_sequences()
        sim.simulate_fiops_v1()

        # Run the evaluation pipeline
        results, fitted_params, elastic_net = pipeline.run_pipeline(
            sim.alignment, sim.attributes_df["fiop"], pipeline.PIPELINE_PARAMS
        )

        params = generalizability_pipeline.GeneralizabilityParams(n_splits=2, n_steps=3)
        (
            model_fit_df,
            model_coef_df,
            feature_counts_df,
        ) = generalizability_pipeline.run_pipeline(
            results, params, pipeline.PIPELINE_PARAMS
        )

        # Plot the results.
        f = generalizability_pipeline.plot_model_fit_vs_library_size(model_fit_df)
        self.assertTrue(f)

        f = generalizability_pipeline.plot_feature_counts_vs_library_size(
            feature_counts_df
        )
        self.assertTrue(f)

        f = generalizability_pipeline.plot_model_hyperparameters_vs_library_size(
            model_fit_df
        )
        self.assertTrue(f)

        f = generalizability_pipeline.plot_feature_counts_vs_model_coefficients(
            model_coef_df, feature_counts_df
        )
        self.assertTrue(f)

        f = generalizability_pipeline.plot_model_coefficients_vs_library_size(
            model_coef_df
        )
        self.assertTrue(f)
