
import pandas as pd
import h2o
from h2o.automl import H2OAutoML
from sklearn.metrics import classification_report, confusion_matrix

# we are loading the data here
train_meta_df = pd.read_csv("GEO_MOUSE_100K_trainin_V2.csv")
test_meta_df  = pd.read_csv("GEO_MOUSE_100K_test_V2.csv")

# the gene expression data
expr_df = pd.read_csv("mouse_expression_data.csv", index_col=0)
expr_df.index.name = "GSM_ID"
expr_df.reset_index(inplace=True)

# selecting target feature
TARGET = "Organ"
train_meta_min = train_meta_df[["GSM_ID", TARGET]]
test_meta_min  = test_meta_df[["GSM_ID", TARGET]]

train_df = pd.merge(train_meta_min, expr_df, on="GSM_ID", how="inner")
test_df  = pd.merge(test_meta_min,  expr_df, on="GSM_ID", how="inner")

print("After merge:")
print("Train shape:", train_df.shape)
print("Test shape:",  test_df.shape)

# removing rows with NA and other unwanted terms 
train_df = train_df.dropna(subset=[TARGET])
test_df  = test_df.dropna(subset=[TARGET])

train_df = train_df[train_df[TARGET] != "Mixed"]
test_df  = test_df[test_df[TARGET] != "Mixed"]

train_df = train_df[train_df[TARGET] != "Normal"]
test_df = test_df[test_df[TARGET] != "Normal"]

print("After dropping rows with missing target:")
print("Train shape:", train_df.shape)
print("Test shape:",  test_df.shape)

exclude_cols = ["GSM_ID", TARGET]
gene_cols = [col for col in train_df.columns if col not in exclude_cols]

train_df = train_df[exclude_cols + gene_cols]
test_df = test_df[exclude_cols + gene_cols]

print("After removing low-count genes:")
print("Train shape:", train_df.shape)
print("Test shape:",  test_df.shape)

#initialize h20
h2o.init()

train_hf = h2o.H2OFrame(train_df)
test_hf  = h2o.H2OFrame(test_df)

all_cols = train_hf.columns
all_cols.remove("GSM_ID")  
all_cols.remove(TARGET)    

features = all_cols
target   = TARGET

# training the model
aml = H2OAutoML(max_runtime_secs=2000, max_models=10, seed=42, stopping_metric="AUTO")
aml.train(x=features, y=target, training_frame=train_hf)

# model evaluation
perf = aml.leader.model_performance(test_hf)
print(perf)

# the code below is used to print model performance
preds = aml.leader.predict(test_hf)
predictions_df = preds.as_data_frame()
true_labels_df = test_hf[TARGET].as_data_frame()
true_labels = true_labels_df[TARGET]
predictions_df['predict'] = predictions_df['predict'].fillna('nan')
true_labels = true_labels.fillna('nan')
predicted_labels = predictions_df['predict']
unique_labels = true_labels.dropna().unique()
report = classification_report(true_labels, predicted_labels, labels=unique_labels, target_names=unique_labels, output_dict=True)
report_df = pd.DataFrame(report).transpose()


print(report_df)

# Save classification report to CSV
report_df.to_csv('classification_report_organ_mousev2.csv')


total_samples = len(true_labels)
correct_predictions = (predicted_labels == true_labels).sum()
accuracy = correct_predictions / total_samples
print(f"Model Accuracy: {accuracy:.4f}")


#output to a file
with open("output_organ_mousev2.txt", "w") as file:
    file.write("Model Performance:\n")
    file.write(str(perf) + "\n\n")

    file.write("Predictions:\n")
    file.write(str(preds.head()) + "\n")
    file.write(f"\nModel Accuracy: {accuracy:.4f}\n")

print("Output written to output_human.txt")

