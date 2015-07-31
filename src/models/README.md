New form of training writes out {model}.{version}.{ext}

feats are the stored features that were used to train the model. 
E.g.
feats.1.0.npz  feats.1.0.pkl  feats.1.0.txt  knn.1.0.pkl  knn.1.0.pkl.stat  knn_cmd.1.0.log  rf.1.0.pkl  rf.1.0.pkl.stat  rf_cmd.1.0.log  svm.1.0.pkl  svm.1.0.pkl.stat  svm_cmd.1.0.log

Notes on Versions:
1.0 - initial NA12878_HG19 without single-end HiC reads
1.1 - NA12878_HG19_GS_0_1 with single-end HiC reads
1.2 - NA12878_HG19_GS_0_1 with single-end HiC reads with SimpleSumFeatures
1.3 - Simulated data on count 001-007, and NA12878_HG_GS_0_1 with SimpleSumFeatures
2.0 - Simulated data on count 001
2.1 - Simulated data on count 004-005
2.2 - Simulated data on count 001-007. 003 and 007 were created on DeBruijn
2.3 - Simulated data on count 001-007 and NA12878_HG_GS_0_1
3.0 - NA12878_HG19_GS_0_1 with single-end HiC reads on WGS-only features
3.1 - NA12878_HG19_GS_0_1 with single-end HiC reads on HiC-only features
3.2 - NA12878_HG19_GS_0_1 and simulated data count 001-007 on WGS-only features
3.3 - NA12878_HG19_GS_0_1 and simulated data count 001-007 on HiC-only features

Repeating Versions to include more extensive correct/predictions per class in outer test stats
1.1 - Same as above, knn, svm, rf
1.3 - Same as above, knn, svm, rf
2.3 - Same as above, knn, svm, rf
3.0 - And later
