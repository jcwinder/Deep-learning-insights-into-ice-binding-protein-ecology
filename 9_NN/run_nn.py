# Setup
import numpy as np
import pandas as pd
from tensorflow.keras import layers
import tensorflow as tf
import matplotlib.pyplot as plt
import sklearn
import keras
import sklearn.datasets
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.utils.class_weight import compute_class_weight
from tensorflow.keras.utils import to_categorical
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras import backend as K
import os
from numpy.random import seed
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
import seaborn as sns

os.environ['TF_XLA_FLAGS'] = '--tf_xla_enable_xla_devices'
seed(132)
tf.random.set_seed(111)

seqfeat_in = pd.read_csv("/path/to/encoded_seqs.csv.gz", compression = 'gzip', low_memory=False) # ESM2-encoded sequences
seqlabels_in = pd.read_csv("/path/to/labeled_seqs.csv") # labelled sequences with a column for names 
seqfeat=pd.DataFrame(seqfeat_in)
seqlabels_df=pd.DataFrame(seqlabels_in)
seqfeat_fl = seqfeat[seqfeat['sequence_id'].isin(seqlabels_df['names'])]
seqlabel_fl = seqlabels_df[seqlabels_df['names'].isin(seqfeat_fl['sequence_id'])]
seqfeat=pd.DataFrame(seqfeat_fl)
print(str(seqfeat.shape))
seqlabel_fl=seqlabel_fl.sort_values('names')
seqfeat = seqfeat.sort_values('sequence_id')
seqlabels_id=pd.DataFrame(seqlabel_fl).loc[:, 'env_id']
seqlabels_id = seqlabels_id.reset_index(drop=True)
seqfeat = seqfeat.reset_index(drop=True).drop(['sequence_id'], axis=1)
seqlbl2=pd.DataFrame(seqlabels_id).loc[:, 'env_id']
seqlbl3=pd.DataFrame(seqlbl2)
target=seqlbl3.iloc[:,:]
features = seqfeat.iloc[:,:]

def dnn_run(n_envs, n_epochs, lrn_rate, batch_size, nlayers, hidden_nodes1, hidden_nodes2=None):
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)

    # Train/validation split
    for initial_train_index, initial_val_index in sss.split(features, target):
        features_train = features.iloc[initial_train_index, :]
        features_val = features.iloc[initial_val_index, :]
        target_train = to_categorical(target.iloc[initial_train_index, :], num_classes=n_envs)
        target_val = to_categorical(target.iloc[initial_val_index, :], num_classes=n_envs)

    # Compute class weights
    y_integers = np.argmax(target_train, axis=1)
    class_weights = compute_class_weight(class_weight='balanced', classes=np.unique(y_integers), y=y_integers)
    d_class_weights = dict(enumerate(class_weights))

    # Reset index for KFold
    features_train = features_train.reset_index(drop=True)
    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    # Initialize lists for storing losses, accuracies, and ROC curve information
    all_train_losses = []
    all_val_losses = []
    all_train_accs = []
    all_val_accs = []
    all_roc_y_true = []
    all_roc_y_pred_probs = []
    all_y_true = []
    all_y_pred = []
    cv_scores = []

    # Loop over each fold
    for fold_train_index, fold_val_index in kf.split(features_train):
        X_train_fold = features_train.iloc[fold_train_index]
        X_val_fold = features_train.iloc[fold_val_index]
        y_train_fold = target_train[fold_train_index]
        y_val_fold = target_train[fold_val_index]
        # Create the model
        model_multi36 = Sequential()
        model_multi36.add(Dense(units=hidden_nodes1, input_dim=X_train_fold.shape[1], activation='relu'))
        
        if nlayers == 2 and hidden_nodes2 is not None:
            model_multi36.add(Dense(units=hidden_nodes2, activation='relu'))

        model_multi36.add(Dense(units=n_envs, activation='softmax'))

        model_multi36.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=lrn_rate),
                              loss='categorical_crossentropy', metrics=['accuracy'])

        # Train the model
        history = model_multi36.fit(X_train_fold, y_train_fold,
                                    epochs=n_epochs,
                                    batch_size=batch_size,
                                    validation_data=(X_val_fold, y_val_fold),
                                    verbose=0,
                                    class_weight=d_class_weights)

        # Store training/validation losses and accuracies
        all_train_losses.append(history.history['loss'])
        all_val_losses.append(history.history['val_loss'])
        all_train_accs.append(history.history['accuracy'])
        all_val_accs.append(history.history['val_accuracy'])

        # Evaluate on the validation set
        y_pred_probs = model_multi36.predict(X_val_fold)
        y_pred_classes = np.argmax(y_pred_probs, axis=1)
        y_true_classes = np.argmax(y_val_fold, axis=1)

        all_y_true.extend(y_true_classes)
        all_y_pred.extend(y_pred_classes)
        cv_scores.append(accuracy_score(y_true_classes, y_pred_classes))


    final_eval = model_multi36.evaluate(x=features_val, y=target_val)
    print("Final evaluation on validation set:", final_eval)

    # Return necessary information for external plotting
    return {
        'train_losses': all_train_losses,
        'val_losses': all_val_losses,
        'train_accs': all_train_accs,
        'val_accs': all_val_accs,
        'initial_train_index': initial_train_index,
        'initial_val_index': initial_val_index,
        'final_eval': final_eval,
        'all_y_true': all_y_true,
        'all_y_pred': all_y_pred,  # Needed for confusion matrix
    }

#Example running
s = dnn_run(n_envs = 5, n_epochs = 100, pca_dimensions = 0, lrn_rate=0.001, batch_size = 256, nlayers = 2, hidden_nodes1 = 800,hidden_nodes2 = 400)

# Assuming `results` is the output of dnn_run
train_losses = s['train_losses']
val_losses = s['val_losses']

# Define specific colors for train and validation curves
train_color = 'blue'
val_color = 'orange'

# Plot learning curves for each fold
for i, (train_loss, val_loss) in enumerate(zip(train_losses, val_losses)):
    # Use the same color for all train curves
    plt.plot(train_loss, label=f'Train Fold {i+1}', color=train_color, linestyle='--')
    # Use the same color for all val curves
    plt.plot(val_loss, label=f'Val Fold {i+1}', color=val_color)

plt.title('Learning Curves')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.show()

class_name_mapping = {
    0: 'frozen_sediment',
    1: 'rock',
    2: 'subsurface',
    3: 'polar_marine',
    4: 'glacier_ice'
}
all_y_true =s['all_y_true']
all_y_pred = s['all_y_pred']

all_y_true_class_names = [class_name_mapping[i] for i in all_y_true]
all_y_pred_class_names = [class_name_mapping[i] for i in all_y_pred]

# Generate confusion matrix
conf_mat = confusion_matrix(all_y_true_class_names, all_y_pred_class_names, labels=sorted(class_name_mapping.values()))

# Plot the confusion matrix using seaborn
plt.figure(figsize=(12, 10))
ax = sns.heatmap(conf_mat / conf_mat.astype(float).sum(axis=1, keepdims=True),
                 fmt='.2%', annot=True, cmap='Blues',
                 xticklabels=sorted(class_name_mapping.values()),
                 yticklabels=sorted(class_name_mapping.values()),
                 vmin=0, vmax=1, annot_kws={'size': 7})
cbar = ax.collections[0].colorbar
cbar.set_ticks([0, .25, 0.5, .75, 1])
cbar.set_ticklabels(['0', '25%', '50%', '75%', '100%'])
plt.xlabel('Predicted')
plt.ylabel('True')
plt.title('Confusion Matrix')
plt.show()









    
    