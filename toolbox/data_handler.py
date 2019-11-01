# Created by woochanghwang at 2019-04-01
"""
list and pandas etc, data handerler
list split
select rows in pandas
"""
import pandas
import pickle
import json
import numpy as np

# Utility function: saves data in JSON format
def dump_json(out_file_name, result):
    with open(out_file_name, 'w') as out_file:
        out_file.write(json.dumps(result, indent=4, separators=(',', ': ')))

# Utility function: loads JSON data into a Python object
def load_json(file_name):
    with open(file_name) as f:
        return json.loads(f.read())

def save_obj(obj, file_addr ):
    with open(file_addr + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(file_addr):
    with open(file_addr+ '.pkl', 'rb') as f:
        return pickle.load(f)

def pandas_select_row_by_list(df,col,selected_list):
    return df[df[col].isin(selected_list)]

def recode_empty_cells(dataframe,fill_value):

    list_of_columns = list(dataframe)
    for column in list_of_columns:
        dataframe = dataframe[column].replace(r'^\s*$', np.nan, regex=True)
        # dataframe[column].fillna(fill_value)

    return dataframe

def pandas_add_column_with_list(df,new,col_name):
    '''

    :param df: dataframe which want to add the col
    :param new: list
    :param col_name: name of col
    :return:
    '''
    # df[col_name] = df.assign(pandas.Series(new).values)
    # return df
    col_name = col_name
    return df.assign(col_name=pandas.Series(new).values)

def pandas_combine_two_df_by_col(df1,df2,col):
    return pandas.merge(df1,df2, on=col)

def list_to_sublist(l,n):
    """
    after receive return
    has to list them
    list(result)
    :param l:
    :param n:
    :return:
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]
def main():
    pass


if __name__ == '__main__':
    main()