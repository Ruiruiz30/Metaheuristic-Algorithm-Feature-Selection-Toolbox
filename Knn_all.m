% K-nearest Neighbor by Jinrui Zhang - (12/31/2024)
% 
function [Acc, MacroFMeasure, MacroPrecision] = Knn_all(feat, label, opts)
% Default of k-value
k = 5;
if isfield(opts,'k'), k = opts.k; end
if isfield(opts,'Model'), Model = opts.Model; end
% Define training & validation sets
trainIdx = Model.training;    testIdx = Model.test;
xtrain   = feat(trainIdx,:);  ytrain  = label(trainIdx);
xvalid   = feat(testIdx,:);   yvalid  = label(testIdx);
% Training model
My_Model = fitcknn(xtrain,ytrain,'NumNeighbors',k); 
% Prediction
pred     = predict(My_Model,xvalid);

% Accuracy
Acc      = sum(pred == yvalid) / length(yvalid);
fprintf('\n Accuracy: %g %%',100 * Acc);

% Calculate per-class precision, recall, and F-measure
numClasses = max(label);
precisionPerClass = zeros(numClasses, 1);
recallPerClass = zeros(numClasses, 1);
FMeasurePerClass = zeros(numClasses, 1);
for c = 1:numClasses
    truePositives = sum((pred == c) & (yvalid == c));
    falsePositives = sum((pred == c) & (yvalid ~= c));
    falseNegatives = sum((pred ~= c) & (yvalid == c));
    precisionPerClass(c) = truePositives / (truePositives + falsePositives);
    recallPerClass(c) = truePositives / (truePositives + falseNegatives);
    FMeasurePerClass(c) = 2 * precisionPerClass(c) * recallPerClass(c) / (precisionPerClass(c) + recallPerClass(c));
end

% Macro-average F-measure
MacroFMeasure = mean(FMeasurePerClass);
MacroPrecision = mean(precisionPerClass);
% Micro-average F-measure
totalTruePositives = sum(truePositives);
totalFalsePositives = sum(falsePositives);
totalFalseNegatives = sum(falseNegatives);
microPrecision = totalTruePositives / (totalTruePositives + totalFalsePositives);
microRecall = totalTruePositives / (totalTruePositives + totalFalseNegatives);
MicroFMeasure = 2 * microPrecision * microRecall / (microPrecision + microRecall);
% F-measure of Macro-average and Micro-average, here we take the macro-average.
% ‌Macro averaging ‌ : Works well on class-imbalanced datasets as it gives equal weight to each class and better reflects performance on small classes.
fprintf('\n Macro-average F-measure: %g', MacroFMeasure);
fprintf('\n Macro-average Precision: %g', MacroPrecision);

end