import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.linear_model import LogisticRegression
from sklearn.utils.class_weight import compute_class_weight
from sklearn.preprocessing import OneHotEncoder
from scipy.sparse import csr_matrix, hstack, vstack
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedShuffleSplit, StratifiedKFold
from sklearn.utils.class_weight import compute_class_weight
from sklearn.metrics import accuracy_score, confusion_matrix

def logreg(X, y):
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
    class_weights = compute_class_weight('balanced', classes=np.unique(y), y=y)
    d_class_weights = dict(zip(np.unique(y), class_weights))

    # Train/validation split
    for initial_train_index, initial_val_index in sss.split(X, y):
        features_train = X[initial_train_index, :]
        features_val = X[initial_val_index, :]
        target_train = y[initial_train_index]
        target_val = y[initial_val_index]

    skf_inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    # Initialize lists for storing accuracies
    logistic_accs = []
    all_y_true = []
    all_y_pred = []
    train_accs = []
    val_accs = []

    target_train_arr = np.asarray(target_train)

    for fold_train_index, fold_val_index in skf_inner.split(features_train, target_train_arr):
        print("Processing fold...")
        X_train_fold = features_train[fold_train_index]
        X_val_fold = features_train[fold_val_index]
        y_train_fold = target_train_arr[fold_train_index]
        y_val_fold = target_train_arr[fold_val_index]
        
        model = LogisticRegression(multi_class="multinomial", solver="saga", max_iter=4000, class_weight=d_class_weights)
        model.fit(X_train_fold, y_train_fold)
    
        # Predict on validation set
        y_train_pred = model.predict(X_train_fold)
        y_val_pred = model.predict(X_val_fold)

        # Store results
        train_accs.append(accuracy_score(y_train_fold, y_train_pred))
        val_accs.append(accuracy_score(y_val_fold, y_val_pred))
        logistic_accs.append(accuracy_score(y_val_fold, y_val_pred))
        all_y_true.extend(y_val_fold)
        all_y_pred.extend(y_val_pred)

    # Plot learning curves
    plt.figure(figsize=(8, 6))
    plt.plot(train_accs, label='Train Accuracy', marker='o')
    plt.plot(val_accs, label='Validation Accuracy', marker='s')
    plt.xlabel('Fold')
    plt.ylabel('Accuracy')
    plt.title('Logistic Regression Learning Curve')
    plt.legend()
    plt.grid()
    plt.savefig("/path/to/logreg_tl.png", format='png')

    # Final evaluation on test set
    print("\nEvaluating on unseen test set...")
    y_test_pred = model.predict(features_val)
    test_acc = accuracy_score(target_val, y_test_pred)
    test_conf_matrix = confusion_matrix(target_val, y_test_pred)
    print(f"Test Set Accuracy: {test_acc:.4f}")
    print("Confusion Matrix (Test Set):\n", test_conf_matrix)

    return {
        'test accuracy': test_acc,
        'test confusion': test_conf_matrix,
        'train accuracies': train_accs,
        'val accuracies': val_accs
    }


seqfeat_in = pd.read_csv("/path/to/sequence_representations_2.csv.gz", compression = 'gzip', low_memory=False) #("/content/sequence_representations_2.csv.gz", compression = "gzip")#("/path/to/sequence_representations.csv", low_memory=False) #"/content/drive/MyDrive/PhD/", low_memory=False)
seqlabels_in = pd.read_csv("/path/to/pf11999_envlbl_7.csv")#("/content/pf11999_envlbl.csv")#("/path/to/pf11999_envlbl.csv")#

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
seqfeat = seqfeat.reset_index(drop=True)
seq_ids = seqfeat['sequence_id'].values  # Save sequence IDs before dropping them
seqfeat = seqfeat.drop(['sequence_id'], axis=1)

seqlbl2=pd.DataFrame(seqlabels_id).loc[:, 'env_id']
seqlbl3=pd.DataFrame(seqlbl2)

target=seqlbl3.iloc[:,:]
features = seqfeat.iloc[:,:]
# Convert features and labels
X = features.to_numpy()
y = target["env_id"]
tst = logreg(X,y)
print("L33 encoded accuracy and confusion matrix")
print(tst['test accuracy'])
print(tst['test confusion'])
print("finished part 1")


OHE data: 
aa_list = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
  # 20 amino acids + gap '-'
# List to store sparse feature matrices for each sequence
alignment_file="/path/to/less215_filtalign.fasta"
encoded_sequences = []
sequence_ids = []

records = SeqIO.parse(alignment_file, "fasta")

for record in SeqIO.parse(alignment_file, "fasta"):
#for record in islice(records, 1):
    #print(record.id)
    sequence_ids.append(record.id)  # Store sequence ID

    # Convert sequence to a DataFrame (each position as a row)
    seq_series = pd.Series(list(str(record.seq)))  

    # One-hot encode using pd.get_dummies()
    oheseq = pd.get_dummies(seq_series, columns=[0], sparse=True).reindex(columns=aa_list, fill_value=0)
    # Convert to a NumPy array and flatten into a single row
    encoded_sequences.append(oheseq.to_numpy().flatten())
    #print(encoded_sequences)
# Convert list of encoded sequences into a DataFrame or NumPy array
X = np.array(encoded_sequences)
# target_df = pd.read_csv("/path/to/pf11999_envlbl_7.csv")
# target_df = target_df[target_df['names'].isin(sequence_ids)]
# target_df = target_df.sort_values('names').reset_index(drop=True)
# y = target_df["env_id"]

np.savetxt('/path/to/ohe_seqs.txt', X, fmt='%d')

X = np.loadtxt('/path/to/ohe_seqs.txt', dtype=int)

alignment_file="/path/to/less215_filtalign.fasta"

sequence_ids = []

records = SeqIO.parse(alignment_file, "fasta")

for record in SeqIO.parse(alignment_file, "fasta"):
#for record in islice(records, 1):
    #print(record.id)
    sequence_ids.append(record.id) 
target_df = pd.read_csv("/path/to/pf11999_envlbl_7.csv")
target_df = target_df[target_df['names'].isin(sequence_ids)]
target_df = target_df.sort_values('names').reset_index(drop=True)
y = target_df["env_id"]

tst2 = logreg(X,y)
print("OHE encoded accuracy and confusion matrix")
print(tst2['test accuracy'])
print(tst2['test confusion'])
print("finished part 2")