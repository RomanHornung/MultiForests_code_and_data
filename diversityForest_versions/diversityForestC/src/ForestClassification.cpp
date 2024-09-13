/*-------------------------------------------------------------------------------
 This file is part of diversityForestC.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of divfor is distributed under MIT license and the
 R package "diversityForestC" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <random>
#include <stdexcept>
#include <cmath>
#include <string>

#include "utility.h"
#include "ForestClassification.h"
#include "TreeClassification.h"
#include "Data.h"

namespace diversityForestC {

void ForestClassification::loadForest(size_t dependent_varID, size_t num_trees,
    std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
    std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values, 
	std::vector<std::vector<size_t>>& forest_split_types, std::vector<std::vector<std::vector<size_t>>>& forest_split_multvarIDs, 
	std::vector<std::vector<std::vector<std::vector<bool>>>>& forest_split_directs, std::vector<std::vector<std::vector<std::vector<double>>>>& forest_split_multvalues,
	std::vector<std::vector<std::vector<size_t>> >& forest_child_muwnodeIDs,
	std::vector<std::vector<std::vector<double>>>& forest_split_muwvalues, 
    std::vector<double>& class_values, std::vector<bool>& is_ordered_variable) {

  this->dependent_varID = dependent_varID;
  this->num_trees = num_trees;
  this->class_values = class_values;
  data->setIsOrderedVariable(is_ordered_variable);

  std::vector<size_t> empty_muw_inds;
  
  // Create trees
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(
        std::make_unique<TreeClassification>(forest_child_nodeIDs[i], forest_split_varIDs[i], forest_split_values[i], forest_split_types[i], forest_split_multvarIDs[i], 
	    forest_split_directs[i], forest_split_multvalues[i], forest_child_muwnodeIDs[i], forest_split_muwvalues[i], empty_muw_inds, &this->class_values, &response_classIDs));
  }

  // Create thread ranges
  equalSplit(thread_ranges, 0, num_trees - 1, num_threads);
}

void ForestClassification::initInternal(std::string status_variable_name) {

  // If npairs not set, use floored square root of number of independent variables.
  if (npairs == 0) {
    unsigned long temp = (size_t)ceil(sqrt((double) (num_variables - 1)) / 2);
    npairs = temp;
  }

  // If mtry not set, use floored square root of number of independent variables.
  if (mtry == 0) {
    unsigned long temp = sqrt((double) (num_variables - 1));
    mtry = std::max((unsigned long) 1, temp);
  }
  
    // If proptry not set, use floored square root of number of independent variables divided by number of independent variables. // asdf
  if (proptry == 0.0) { // asdf
    double temp = sqrt((double) (num_variables - 1)) / (double) (num_variables - 1); // asdf
    proptry = std::min((double) 1, temp); // asdf
  } // asdf

  // Set minimal node size
  if (min_node_size == 0) {
    min_node_size = DEFAULT_MIN_NODE_SIZE_CLASSIFICATION;
  }

  // Create class_values and response_classIDs
  if (!prediction_mode) {
    for (size_t i = 0; i < num_samples; ++i) {
      double value = data->get(i, dependent_varID);

      // If classID is already in class_values, use ID. Else create a new one.
      uint classID = find(class_values.begin(), class_values.end(), value) - class_values.begin();
      if (classID == class_values.size()) {
        class_values.push_back(value);
      }
      response_classIDs.push_back(classID);
    }
  }

  // Create sampleIDs_per_class if required
  if (sample_fraction.size() > 1) {
    sampleIDs_per_class.resize(sample_fraction.size());
    for (auto& v : sampleIDs_per_class) {
      v.reserve(num_samples);
    }
    for (size_t i = 0; i < num_samples; ++i) {
      size_t classID = response_classIDs[i];
      sampleIDs_per_class[classID].push_back(i);
    }
  }

  // Set class weights all to 1
  class_weights = std::vector<double>(class_values.size(), 1.0);

  // Sort data if memory saving mode
  if (!memory_saving_splitting) {
    data->sort();
  }
}

void ForestClassification::growInternal() {
  trees.reserve(num_trees);
  for (size_t i = 0; i < num_trees; ++i) {
    trees.push_back(
        std::make_unique<TreeClassification>(&class_values, &response_classIDs, &sampleIDs_per_class, &class_weights));
  }
}

void ForestClassification::allocatePredictMemory() {
  size_t num_prediction_samples = data->getNumRows();
  if (predict_all || prediction_type == TERMINALNODES) {
    predictions = std::vector<std::vector<std::vector<double>>>(1,
        std::vector<std::vector<double>>(num_prediction_samples, std::vector<double>(num_trees)));
  } else {
    predictions = std::vector<std::vector<std::vector<double>>>(1,
        std::vector<std::vector<double>>(1, std::vector<double>(num_prediction_samples)));
  }
}

void ForestClassification::predictInternal(size_t sample_idx) {
  if (predict_all || prediction_type == TERMINALNODES) {
    // Get all tree predictions
    for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
      if (prediction_type == TERMINALNODES) {
        predictions[0][sample_idx][tree_idx] = getTreePredictionTerminalNodeID(tree_idx, sample_idx);
      } else {
		if (divfortype == 1) {
          predictions[0][sample_idx][tree_idx] = getTreePrediction(tree_idx, sample_idx);
		}
		if (divfortype == 2) {
          predictions[0][sample_idx][tree_idx] = getTreePredictionMultivariate(tree_idx, sample_idx);
		}
		if (divfortype == 3) {
          predictions[0][sample_idx][tree_idx] = getTreePredictionMuw(tree_idx, sample_idx);
		}
      }
    }
  } else {
    // Count classes over trees and save class with maximum count
    std::unordered_map<double, size_t> class_count;
    for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
	  if (divfortype == 1) {
        ++class_count[getTreePrediction(tree_idx, sample_idx)];
	  }
	  if (divfortype == 2) {
		++class_count[getTreePredictionMultivariate(tree_idx, sample_idx)];
	  }
	  if (divfortype == 3) {
		++class_count[getTreePredictionMuw(tree_idx, sample_idx)];
	  }
    }
    predictions[0][0][sample_idx] = mostFrequentValue(class_count, random_number_generator);
  }
}

void ForestClassification::computePredictionErrorInternal() {

  // Class counts for samples
  std::vector<std::unordered_map<double, size_t>> class_counts;
  class_counts.reserve(num_samples);
  for (size_t i = 0; i < num_samples; ++i) {
    class_counts.push_back(std::unordered_map<double, size_t>());
  }

  // For each tree loop over OOB samples and count classes
  for (size_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    for (size_t sample_idx = 0; sample_idx < trees[tree_idx]->getNumSamplesOob(); ++sample_idx) {
	  size_t sampleID = trees[tree_idx]->getOobSampleIDs()[sample_idx];
	  if (divfortype == 1) {
        ++class_counts[sampleID][getTreePrediction(tree_idx, sample_idx)];
	  }
	  if (divfortype == 2) {
        ++class_counts[sampleID][getTreePredictionMultivariate(tree_idx, sample_idx)];
	  }
	  if (divfortype == 3) {
        ++class_counts[sampleID][getTreePredictionMuw(tree_idx, sample_idx)];
	  }
    }
  }

  // Compute majority vote for each sample
  predictions = std::vector<std::vector<std::vector<double>>>(1,
      std::vector<std::vector<double>>(1, std::vector<double>(num_samples)));
  for (size_t i = 0; i < num_samples; ++i) {
    if (!class_counts[i].empty()) {
      predictions[0][0][i] = mostFrequentValue(class_counts[i], random_number_generator);
    } else {
      predictions[0][0][i] = NAN;
    }
  }

  // Compare predictions with true data
  size_t num_missclassifications = 0;
  size_t num_predictions = 0;
  for (size_t i = 0; i < predictions[0][0].size(); ++i) {
    double predicted_value = predictions[0][0][i];
    if (!std::isnan(predicted_value)) {
      ++num_predictions;
      double real_value = data->get(i, dependent_varID);
      if (predicted_value != real_value) {
        ++num_missclassifications;
      }
      ++classification_table[std::make_pair(real_value, predicted_value)];
    }
  }
  overall_prediction_error = (double) num_missclassifications / (double) num_predictions;
}

// #nocov start
void ForestClassification::saveToFileInternal(std::ofstream& outfile) {

  // Write num_variables
  outfile.write((char*) &num_variables, sizeof(num_variables));

  // Write treetype
  TreeType treetype = TREE_CLASSIFICATION;
  outfile.write((char*) &treetype, sizeof(treetype));

  // Write class_values
  saveVector1D(class_values, outfile);
}

double ForestClassification::getTreePrediction(size_t tree_idx, size_t sample_idx) const {
  const auto& tree = dynamic_cast<const TreeClassification&>(*trees[tree_idx]);
  return tree.getPrediction(sample_idx);
}

double ForestClassification::getTreePredictionMultivariate(size_t tree_idx, size_t sample_idx) const {
  const auto& tree = dynamic_cast<const TreeClassification&>(*trees[tree_idx]);
  return tree.getPredictionMultivariate(sample_idx);
}

double ForestClassification::getTreePredictionMuw(size_t tree_idx, size_t sample_idx) const {
  const auto& tree = dynamic_cast<const TreeClassification&>(*trees[tree_idx]);
  return tree.getPredictionMuw(sample_idx);
}

size_t ForestClassification::getTreePredictionTerminalNodeID(size_t tree_idx, size_t sample_idx) const {
  const auto& tree = dynamic_cast<const TreeClassification&>(*trees[tree_idx]);
  return tree.getPredictionTerminalNodeID(sample_idx);
}

// Multi Forests: Compute variable importance.
void ForestClassification::computeImportanceMuw() {

  // Compute EIM values in multiple threads
  progress = 0;
#ifdef R_BUILD
  aborted = false;
  aborted_threads = 0;
#endif

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  // Initialize importance
  std::vector<std::vector<double>> var_imp_multiway_threads(num_threads);
  std::vector<std::vector<double>> var_imp_binary_threads(num_threads);

  // Compute importance
  for (uint i = 0; i < num_threads; ++i)
  {    
    if (importance_mode == MUWIMP_MULTIWAY || importance_mode == MUWIMP_BOTH)
    {
       var_imp_multiway_threads[i].resize(metricind.size(), 0);
    }
    if (importance_mode == MUWIMP_BINARY || importance_mode == MUWIMP_BOTH)
    {
       var_imp_binary_threads[i].resize(num_independent_variables, 0);
    }
    threads.emplace_back(&ForestClassification::computeTreeImportanceMuwInThread, this, i,
                         std::ref(var_imp_multiway_threads[i]), std::ref(var_imp_binary_threads[i]));
  }
  showProgress("Computing variable importance..", num_trees);
  for (auto &thread : threads)
  {
    thread.join();
  }
  
#ifdef R_BUILD
  if (aborted_threads > 0)
  {
    throw std::runtime_error("User interrupt.");
  }
#endif

  // Sum thread importances
  if (importance_mode == MUWIMP_MULTIWAY || importance_mode == MUWIMP_BOTH)
  {
	  var_imp_multiway.resize(metricind.size(), 0);
  for (size_t i = 0; i < metricind.size(); ++i)
  {
    for (uint j = 0; j < num_threads; ++j)
    {
      var_imp_multiway[i] += var_imp_multiway_threads[j][i];
    }
  }  
  }
  
    if (importance_mode == MUWIMP_BINARY || importance_mode == MUWIMP_BOTH)
  {
	  var_imp_binary.resize(num_independent_variables, 0);
  for (size_t i = 0; i < num_independent_variables; ++i)
  {
    for (uint j = 0; j < num_threads; ++j)
    {
      var_imp_binary[i] += var_imp_binary_threads[j][i];
    }
  }  
  }

  var_imp_multiway_threads.clear();
  var_imp_binary_threads.clear();

  if (importance_mode == MUWIMP_MULTIWAY || importance_mode == MUWIMP_BOTH)
  {
  for (size_t i = 0; i < var_imp_multiway.size(); ++i)
  {
    var_imp_multiway[i] /= num_trees;
  }
  }

  if (importance_mode == MUWIMP_BINARY || importance_mode == MUWIMP_BOTH)
  {
  for (size_t i = 0; i < var_imp_binary.size(); ++i)
  {
    var_imp_binary[i] /= num_trees;
  }
  }

}

// Multi Forests: Compute variable importance in one thread.
void ForestClassification::computeTreeImportanceMuwInThread(uint thread_idx, std::vector<double>& importance_multiway,
     std::vector<double>& importance_binary) {
  if (thread_ranges.size() > thread_idx + 1) {
    for (size_t i = thread_ranges[thread_idx]; i < thread_ranges[thread_idx + 1]; ++i) {
      if (auto treeClassPtr = dynamic_cast<TreeClassification*>(trees[i].get())) {
        treeClassPtr->computeImportanceMuw(importance_multiway, importance_binary);
      }

      // Check for user interrupt
#ifdef R_BUILD
      if (aborted) {
        std::unique_lock<std::mutex> lock(mutex);
        ++aborted_threads;
        condition_variable.notify_one();
        return;
      }
#endif

      // Increase progress by 1 tree
      std::unique_lock<std::mutex> lock(mutex);
      ++progress;
      condition_variable.notify_one();
    }
  }
}

// #nocov end

}// namespace diversityForestC
