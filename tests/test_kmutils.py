import pytest
import lifelines
import pandas as pd

from io import StringIO
from vtg_km.v6_km_utils import get_km_event_table, get_unique_event_times


class TestKMUtils:
    @pytest.fixture(scope="function")
    def lifeline_waltons(self):
        waltons = lifelines.datasets.load_waltons()
        # Cheaply and badly check if the dataset has changed, as we want to
        # further make sure differences in test results are due to potential
        # errors in our code, not unexpected upstream changes to datasets from
        # lifeines.
        if pd.util.hash_pandas_object(waltons).sum() != 11603055737657237860:
            pytest.exit("Waltons dataset from lifelines _seems_ to have changed. Aborting tests with it.")

        kmf = lifelines.KaplanMeierFitter()
        kmf.fit(waltons['T'], waltons['E'])

        return waltons, kmf

    def test_get_km_event_table_waltons_no_bining(self, lifeline_waltons):
        waltons_df, kmf_waltons = lifeline_waltons
        params = {
            "time_column_name": "T",
            "unique_event_times": waltons_df['T'].unique().tolist(),
            "censor_column_name": "E",
            "bin_size": None,
        }

        our_table_json = get_km_event_table(mock_data=[waltons_df], **params)
        our_table = pd.read_json(StringIO(our_table_json))
        lifelines_table = kmf_waltons.event_table

        # lifelines' uses event_at for index, with no extra index, we adapt ours
        our_table.set_index('T', inplace=True)

        # lifelines' has some extra columns we don't generate
        lifelines_table.drop(columns=["entrance"], inplace=True)

        # TODO: Check if this is something we should do
        # lifelines' adds a row for time 0, which we don't
        lifelines_table.drop(index=0, inplace=True)

        # lifelines' uses float, since integers are a proper subset of floats,
        # we convert ours
        our_table.index = pd.to_numeric(our_table.index, errors='coerce').astype(float)

        # Check if the two event tables are the same
        assert our_table.equals(lifelines_table)

    @pytest.mark.skip(reason="Not implemented yet")
    def test_get_km_event_table_waltons_with_binning(self, lifeline_waltons):
        pass

    def test_get_unique_event_times(self, lifeline_waltons):
        waltons_df, kmf_waltons = lifeline_waltons
        params = {
            "time_column_name": "T"
        }
        our_unique_times = sorted(get_unique_event_times(mock_data=[waltons_df], **params))
        lifelines_times = kmf_waltons.event_table.drop(index=0).index.to_list()

        assert our_unique_times == lifelines_times

