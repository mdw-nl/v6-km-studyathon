import time
import numpy as np
import pandas as pd
from vantage6.tools.util import info


def master(
        client, data, time_col, censor_col, method='km', bins=None,
        organization_ids=None
):
    """This package does the following:
            2. Calculates the coordinates of the Kaplan Meier curve
    """

    info('Collecting participating organizations')
    if isinstance(organization_ids, list) is False:
        organizations = client.get_organizations_in_my_collaboration()
        ids = [organization.get("id") for organization in organizations]
    else:
        ids = organization_ids

    info(f'Sending task to organizations {ids}')
    km, local_event_tables = calculate_KM(
        client, ids, time_col, censor_col, method, bins
    )

    return {'kaplanMeier': km, 'local_event_tables': local_event_tables}


def calculate_KM(client, ids, time_col, censor_col, method, bins):

    if method == 'km':
        # Get local unique time events
        kwargs_dict = {'time_col': time_col}
        method = 'get_unique_event_times'
        results = subtaskLauncher(client, [method, kwargs_dict, ids])

        # Combine local unique times into global times
        local_uet = []
        for output in results:
            local_uet.append(output['unique_event_times'])
        local_uet.append([0])
        unique_event_times = list(
            set([item for sublist in local_uet for item in sublist])
        )
    elif method == 'binning':
        try:
            # Define bins for time events
            unique_event_times = list(
                range(0, bins['maxtime']+bins['size'], bins['size'])
            )
        except Exception as e:
            info(f'Exception occurred with input \'bins\': {e}')
    else:
        info(f'Unknown method: {method}')

    ##### 2) Ask to calculate local event tables #####

    kwargs_dict = {
        'time_col': time_col,
        'censor_col': censor_col,
        'unique_event_times': unique_event_times,
        'method': method
    }
    method = 'get_km_event_table'
    results = subtaskLauncher(client, [method, kwargs_dict, ids])

    local_event_tables = []
    for output in results:
        local_event_tables.append(output['event_table'])

    km = pd.concat(local_event_tables).groupby(time_col, as_index=False).sum()
    km['Hazard'] = km['Deaths'] / km['AtRisk']
    km['Surv'] = (1 - km['Hazard']).cumprod()
    km['cdf'] = 1 - km['Surv']
    km['pmf'] = np.diff(km['cdf'], prepend=0)
    return km, local_event_tables




def RPC_get_unique_event_times(
        data, time_col, data_set=None, filt=None, median_lp=None
):
    """Get Unique Event Times
    """
    # df = data_selector(data, data_set,filt, median_lp) #data_set, filt=None, median_lp=None)
    return {
        "unique_event_times": data[time_col].unique()
    }


def RPC_get_km_event_table(
        data, time_col, unique_event_times, censor_col, method
):
    df = data.copy()

    if method == 'binning':
        # Bin event time data
        df[time_col] = np.float16(pd.cut(
            df[time_col], bins=unique_event_times, labels=unique_event_times[1:]
        ))

    info(str(len(df)))
    death = df.groupby(time_col, as_index=False).sum().rename(columns={censor_col: 'Deaths'})[[time_col, 'Deaths']]
    death = pd.DataFrame(unique_event_times, columns=[time_col]).merge(death, on=time_col, how='left').fillna(0)

    total = df.groupby(time_col, as_index=False).count().rename(columns={censor_col: 'Total'})[[time_col, 'Total']]
    total = pd.DataFrame(unique_event_times, columns=[time_col]).merge(total, on=time_col, how='left').fillna(0)

    km = death.merge(total, on=time_col)

    km = km.merge(pd.DataFrame.from_dict({unique_event_times[i]: len(df[df[time_col] >= unique_event_times[i]]) for
                                          i in range(len(unique_event_times))}, orient='index').rename(
        columns={0: 'AtRisk'}).sort_index(), left_on=time_col, right_index=True)
    return {'event_table':km}


def subtaskLauncher(client, taskInfo):
    method, kwargs_dict, ids = taskInfo

    info(f'sending task to organizations {ids}')

    task = client.create_new_task(
        input_={
            'method': method,
            'kwargs': kwargs_dict
        },
        organization_ids=ids
    )

    info("Waiting for results")
    task_id = task.get("id")
    task = client.get_task(task_id)
    while not task.get("complete"):
        task = client.get_task(task_id)
        info("Waiting for results")
        time.sleep(1)
    # Once we know the partials are complete, we can collect them.
    info("Obtaining results")
    results = client.get_results(task_id=task.get("id"))
    return results
