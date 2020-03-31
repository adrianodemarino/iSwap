
import pandas as pd 
import numpy as np 
import random as rd
import math 
import matplotlib.pyplot as plt
import scipy.special as sps

hg = pd.read_csv("hg38.size",sep="\t")
sample = pd.read_csv("sample_name.tsv",sep="\t")

#Sequence of Interest Sofi Creator
def Sofi():
	selection = hg.sample(1)
	chrom = selection.CHR.values
	position = rd.randrange(selection.Size.values)
	Construct = chrom+"_"+str(position)
	return Construct[0]

def GammaDistribution():	
	shape, scale = 2., 2. 
	s = np.random.gamma(shape, scale, 1000)
	count, bins, ignored = plt.hist(s, 50, density=True)
	y = bins**(shape-1)*(np.exp(-bins/scale)/(sps.gamma(shape)*scale**shape))
	plt.plot(bins, y, linewidth=2, color='r')
	plt.show()

def ShearsiteDistribution(len=1000000):
	shape, scale = 2., 2. 
	s = np.random.gamma(shape, scale, len)
	s = (s*100).astype(int)
	s = s[s<=1000]
	return s

order_i = {1:10,2:100,3:1000,4:10000,5:100000,6:1000000,7:10000000,8:100000000,9:1000000000,10:10000000000}
order_f = {1:99,2:999,3:9999,4:99999,5:999999,6:9999999,7:99999999,8:999999999,9:999999999,10:9999999999} 

def magnitude (value): 
    if (value == 0): return 0 
    return int(math.floor(math.log10(abs(value))))

def getSimuNumber (value): 
    real_magnitude = magnitude(value) 
    swap_magnitude = real_magnitude - 2
    return np.random.randint(1*order_i[swap_magnitude],1*order_f[swap_magnitude]) 

df["swap_magnitude"] = df.swapped_reads.apply(lambda x: magnitude(x))
df["reads_magnitude"] = df["#reads"].apply(lambda x: magnitude(x))

df['barcode'] = df['barcode'].map({1: 'HIGH', 2: 'MEDIUM', 3: 'LOW'}) 
df = pd.get_dummies(df, prefix='', prefix_sep='') 

df = pd.read_csv("to_predict.tsv",sep="\t")
train_df = df.sample(frac=0.8,random_state=0)
test_df = df.drop(train_df.index)
#_ =sns.pairplot(train_df[["swapped_reads", "#reads","swap_magnitude","reads_magnitude","barcode"]], diag_kind="kde")
train_stats = train_df.describe()
train_stats.pop("swapped_reads")
train_stats = train_stats.transpose()
train_labels = train_df.pop('swapped_reads')
test_labels = test_df.pop('swapped_reads')
def norm(x):
  return (x - train_stats['mean']) / train_stats['std']
normed_train_data = norm(train_df)
normed_test_data = norm(test_df)
normed_test_data
def build_model():
  model = keras.Sequential([
    layers.Dense(64, activation='relu', input_shape=[len(train_df.keys())]),
    layers.Dense(64, activation='relu'),
    layers.Dense(1)
  ])

  optimizer = tf.keras.optimizers.RMSprop(0.001)

  model.compile(loss='mse',
                optimizer=optimizer,
                metrics=['mae', 'mse'])
  return model
model = build_model()
model.summary()
example_batch = normed_train_data[:10]
example_result = model.predict(example_batch)
example_result
EPOCHS = 1000

history = model.fit(
  normed_train_data, train_labels,
  epochs=EPOCHS, validation_split = 0.2, verbose=0,
  callbacks=[tfdocs.modeling.EpochDots()])
hist = pd.DataFrame(history.history)
hist['epoch'] = history.epoch
hist.tail()
plotter = tfdocs.plots.HistoryPlotter(smoothing_std=2)
plotter.plot({'Basic': history}, metric = "mae")
plotter.plot({'Basic': history}, metric = "mse")
model = build_model()

# The patience parameter is the amount of epochs to check for improvement
early_stop = keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)

early_history = model.fit(normed_train_data, train_labels, 
                    epochs=EPOCHS, validation_split = 0.2, verbose=0, 
                    callbacks=[early_stop, tfdocs.modeling.EpochDots()])
plotter.plot({'Early Stopping': early_history}, metric = "mae")
plotter.plot({'Early Stopping': early_history}, metric = "mae")
loss, mae, mse = model.evaluate(normed_test_data, test_labels, verbose=2)
print("Testing set Mean Abs Error: {:5.2f} swapped_reads".format(mae))
test_predictions = model.predict(normed_test_data).flatten()
a = plt.axes(aspect='equal')
_ = plt.scatter(test_labels, test_predictions)
error = test_predictions - test_labels
_ = plt.hist(error, bins = 25)
