In this study, we built trinary classification models to predict an
organism’s genetic composition based on their codon usage
frequencies. The three classes were nuclear, mitochondrial and
chloroplast DNA. Misclassification error costs, the costs of making
various types of classification errors, were assumed to be equal. Two
sets of models were built, one with the original training data and
another with training data that was transformed (square root
transformation). Classification and regression trees, tree-bag, random
forests, C5.0, k-nearest neighbor, and neural networks models were
built.
Order of using files:
1) codon_0_data_imputation.R
2) codon_1_eda.R
3) codon_2_train_test_split.R
4) codon_3a_train_trans.R
5) codon_3a_model_training_NOT_processed.R
6) codon_3a_model_training_processed.R
