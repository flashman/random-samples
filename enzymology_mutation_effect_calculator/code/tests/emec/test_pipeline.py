import json
import pathlib
import tempfile
import unittest

from emec.pipeline import run_pipeline, PIPELINE_PARAMS
from emec.utils.simulation_utils import RandomMutagenesisSimulation
from emec.utils.statsmodels_utils import get_model_coefficients

from session import MutationEffectCalculatorSession


import tests.fixtures.fern_gar0722_ssl_small_2021_08_04 as ssl_small


class TestPipeline(unittest.TestCase):
    def test_run_pipline_degree_one(self):
        # For this test, ensure that beneficial mutations are always well separated from neutral
        # mutations.

        mutation_params = {"n_mutations": 3}
        fiop_params = {"n_beneficial": 5, "benefit_mu": 2.0, "benefit_sigma": 0.05}

        sim = RandomMutagenesisSimulation(10, 2000, mutation_params, fiop_params)
        sim.simulate_sequences()
        sim.simulate_fiops_v1()

        # Run the evaluation pipeline
        results, fitted_params, elastic_net = run_pipeline(
            sim.alignment, sim.attributes_df["fiop"], PIPELINE_PARAMS
        )

        # Extract predicted beneficial mutation names from the fitted model.
        # Only choose mutations that have "large" coefficients.
        results_df = get_model_coefficients(results)
        found_muts = set(results_df[results_df.coef > 0.5].index) - {"const"}

        # Extract beneficial mutations from the simulation.
        expected_muts = set(sim.model_df.name.unique())

        # Check that all the expected beneficial mutations were found..
        self.assertSetEqual(expected_muts, found_muts)

    def test_run_pipline_degree_two(self):
        # For this test, ensure that beneficial mutations are always well separated from neutral
        # mutations.

        mutation_params = {"n_mutations": 2}
        fiop_params = {"n_beneficial": 5, "benefit_mu": 2.0, "benefit_sigma": 0.05}

        sim = RandomMutagenesisSimulation(10, 1000, mutation_params, fiop_params)
        sim.simulate_sequences()
        sim.simulate_fiops_v1()

        # Run the evaluation pipeline
        pipeline_params = PIPELINE_PARAMS.copy()
        pipeline_params["degree"] = 2
        pipeline_params["elastic_net_cv"]["l1_ratio"] = [0.1, 1.0]
        pipeline_params["group_k_fold"]: {"n_splits": 3}

        results, fitted_params, elastic_net = run_pipeline(
            sim.alignment, sim.attributes_df["fiop"], pipeline_params
        )

        # Extract predicted beneficial mutation names from the fitted model.
        # Only choose mutations with significant coefficients.
        # While strictly wrong allow interaction terms for now.
        results_df = get_model_coefficients(results)
        found_muts = set(results_df[results_df["P>|t|"] < 0.00001].index) - {"const"}
        found_muts = set(m for mm in found_muts for m in mm.split())

        # Extract beneficial mutations from the simulation.
        expected_muts = set(sim.model_df.name.unique())

        # Check that all the expected beneficial mutations were found..
        self.assertSetEqual(expected_muts, found_muts)

    def test_run_pipline_ssl(self):
        # Use session object to simplify fixtures imports.
        # TODO: Add SSL data simulation.
        with tempfile.TemporaryDirectory() as output_dir:
            session = MutationEffectCalculatorSession(
                input_dir=ssl_small.FIXTURES_PATH, output_dir=output_dir
            )
            session.load(
                performance_filename=pathlib.Path(ssl_small.ATTRIBUTES_PATH).name,
                seq_filename=pathlib.Path(ssl_small.DNA_SEQUENCE_PATH).name,
                reference_seq_id=None,
                feature_term=None,
                annotation_name=None,
            )

            # Run the pipeline.
            results, fitted_params, elastic_net = run_pipeline(
                session.alignment,
                session.attributes_df[ssl_small.FIOP_COLUMNE],
                {"precondition": True},
            )

            # Extract predicted beneficial mutation names from the fitted model.
            coef_df = get_model_coefficients(results, p_threshold=0.000001)
            beneficial = list(coef_df[coef_df.coef > 0].index)
            deleterious = list(coef_df[coef_df.coef < 0].index)

            # There is no ground truth for this dataset.  So just use the experimental results.
            expected_beneficial = ["const", "D103E", "Y153E", "F220G"]
            expected_deletarious = [
                "A86E",
                "I91K",
                "A132F",
                "I149E",
                "L181G",
                "L211E",
                "E218F",
                "L253E",
            ]

            assert beneficial == expected_beneficial
            assert deleterious == expected_deletarious
