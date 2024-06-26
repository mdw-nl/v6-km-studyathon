import warnings
import pandas as pd
from vantage6.algorithm.tools.mock_client import MockAlgorithmClient
warnings.filterwarnings("ignore")


# Initialize the mock server. The datasets simulate the local datasets from
# the node. In this case we have two parties having two different datasets:
# a.csv and b.csv. The module name needs to be the name of your algorithm
# package. This is the name you specified in `setup.py`, in our case that
# would be v6-correlation-matrix-py.
dataset_1 = {"database": "./vtg_km/local/data1.csv", "db_type": "csv"}
dataset_2 = {"database": "./vtg_km/local/data2.csv", "db_type": "csv"}
org_ids = ids = [0, 1]

client = MockAlgorithmClient(
    datasets = [[dataset_1], [dataset_2]],
    organization_ids=org_ids,
    module="vtg_km"
)

organizations = client.organization.list()
org_ids = ids = [organization["id"] for organization in organizations]

# To trigger the master method you also need to supply the `master`-flag
# to the input. Also note that we only supply the task to a single organization
# as we only want to execute the central part of the algorithm once. The master
# task takes care of the distribution to the other parties.
average_task = client.task.create(
    input_={
        'method': 'master',
        'kwargs': {
            'time_column_name': 'TIME_AT_RISK',
            'censor_column_name': 'MORTALITY_FLAG',
            'bin_size': 5,
            'query_string': 'COHORT_DEFINITION_ID == 1029'
        }
    },
    organizations=[org_ids[0]]
)

results = client.result.get(average_task.get("id"))
df_events = pd.read_json(results)
print(df_events.head())
