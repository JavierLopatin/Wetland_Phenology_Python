
import pandas as pd

df = pd.read_csv('data/clusters_rmse.csv')
df.head()

# get 5000 random rows from the dataframe
df_random = df.sample(n=5000, random_state=1)

df_random.to_csv('data/rmse_random.csv', index=False)