/*-------------------------------------------------------------------------------
 This file is part of diversityForestD.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of divfor is distributed under MIT license and the
 R package "diversityForestD" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#ifndef TREEPROBABILITY_H_
#define TREEPROBABILITY_H_

#include <map>
#include <vector>

#include "globals.h"
#include "Tree.h"

namespace diversityForestD {

class TreeProbability: public Tree {
public:
  TreeProbability(std::vector<double>* class_values, std::vector<uint>* response_classIDs,
      std::vector<std::vector<size_t>>* sampleIDs_per_class, std::vector<double>* class_weights);

  // Create from loaded forest
  TreeProbability(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
      std::vector<double>& split_values, std::vector<size_t>& split_types, std::vector<std::vector<size_t>>& split_multvarIDs, 
	  std::vector<std::vector<std::vector<bool>>>& split_directs, std::vector<std::vector<std::vector<double>>>& split_multvalues, std::vector<std::vector<size_t>>& child_muwnodeIDs, 
	  std::vector<std::vector<double>>& split_muwvalues, std::vector<size_t>& muw_inds, std::vector<double>* class_values, std::vector<uint>* response_classIDs,
      std::vector<std::vector<double>>& terminal_class_counts);

  TreeProbability(const TreeProbability&) = delete;
  TreeProbability& operator=(const TreeProbability&) = delete;

  virtual ~TreeProbability() override = default;

  void allocateMemory() override;

  void grow(std::vector<double>* variable_importance) override;

  void addToTerminalNodes(size_t nodeID);
  void computePermutationImportanceInternal(std::vector<std::vector<size_t>>* permutations);
  void appendToFileInternal(std::ofstream& file) override;

  void predictMuw(const Data* prediction_data, bool oob_prediction) override;

  void computeImportanceMuw(std::vector<double> &forest_multiway, std::vector<double> &forest_binary);
  void dropDownRecordFirstVisited(std::vector<std::vector<int>> &all_OOBs_first_varIDs, size_t muw_ind);
  double computeImportanceDifference(std::vector<size_t> oob_sampleIDs_subset, std::vector<int> first_visited_nodeID_oob, size_t muw_ind);
  double computeImportanceNode(size_t nodeID, std::vector<size_t> oob_sampleIDs_subset_subset, size_t muw_ind);
  double computeImportanceNodePermuted(size_t nodeID, std::vector<size_t> oob_sampleIDs_subset_subset, std::vector<size_t> oob_sampleIDs_subset_subset_permuted, size_t muw_ind);
  
  const std::vector<double>& getPrediction(size_t sampleID) const {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[sampleID];
    return terminal_class_counts[terminal_nodeID];
  }

  size_t getPredictionTerminalNodeID(size_t sampleID) const {
    return prediction_terminal_nodeIDs[sampleID];
  }

  const std::vector<std::vector<double>>& getTerminalClassCounts() const {
    return terminal_class_counts;
  }
  
  const std::vector<std::vector<double>>& getSplitMuwValues() const {
    return split_muwvalues;
  }
    
  const std::vector<std::vector<size_t>>& getChildMuwNodeIDs() const {
    return child_muwnodeIDs;
  }
  
  const std::vector<size_t>& getMuwInds() const {
    return muw_inds;
  }  

private:
  bool splitNode(size_t nodeID) override;
  bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) override;
  bool splitNodeUnivariateInternal(size_t nodeID, std::vector<std::pair<size_t, double>> sampled_varIDs_values) override; // asdf
  bool checkWhetherFinal(size_t nodeID, std::vector<size_t>& varIDs_rel);
  bool splitNodeMultivariateInternal(size_t nodeID, std::vector<size_t> sampled_split_types, std::vector<std::vector<size_t>> sampled_split_multvarIDs, std::vector<std::vector<std::vector<bool>>> sampled_split_directs, std::vector<std::vector<std::vector<double>>> sampled_split_multvalues) override;
  void drawSplitsMuw(size_t nodeID, std::vector<std::vector<double>>& split_muwvalues_temp, std::vector<size_t>& varIDs_temp, std::vector<size_t> varIDs_rel);
  void splitNodeMuwMuwInternal(size_t nodeID, std::vector<std::vector<double>> split_muwvalues_temp, std::vector<size_t> varIDs_temp);
  void splitNodeMuwUnivInternal(size_t nodeID, std::vector<std::vector<double>> split_muwvalues_temp, std::vector<size_t> varIDs_temp);
  void createEmptyNodeInternal() override;

  double computePredictionAccuracyInternal() override;

  // Called by splitNodeInternal(). Sets split_varIDs and split_values.
  bool findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  bool findBestSplitUnivariate(size_t nodeID, std::vector<std::pair<size_t, double>> sampled_varIDs_values); // asdf
  bool findBestSplitMultivariate(size_t nodeID, std::vector<size_t> sampled_split_types, std::vector<std::vector<size_t>> sampled_split_multvarIDs, std::vector<std::vector<std::vector<bool>>> sampled_split_directs, std::vector<std::vector<std::vector<double>>> sampled_split_multvalues);
  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease);
  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease, const std::vector<double>& possible_split_values, std::vector<size_t>& class_counts_right,
      std::vector<size_t>& n_right);
  void findBestSplitValueLargeQ(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease);
  void findBestSplitValueUnordered(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease);

  bool findBestSplitExtraTrees(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void findBestSplitValueExtraTrees(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease);
  void findBestSplitValueExtraTrees(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease, const std::vector<double>& possible_split_values, std::vector<size_t>& class_counts_right,
      std::vector<size_t>& n_right);
  void findBestSplitValueExtraTreesUnordered(size_t nodeID, size_t varID, size_t num_classes,
      const std::vector<size_t>& class_counts, size_t num_samples_node, double& best_value, size_t& best_varID,
      double& best_decrease);

  void addImpurityImportance(size_t nodeID, size_t varID, double decrease);

  void bootstrapClassWise() override;
  void bootstrapWithoutReplacementClassWise() override;

  void cleanUpInternal() override {
    counter.clear();
    counter.shrink_to_fit();
    counter_per_class.clear();
    counter_per_class.shrink_to_fit();
  }

  // Multi forests: Vector of child node IDs (the i-th element contains
  // the indices of the nodes that are child nodes of node i):
  std::vector<std::vector<size_t>> child_muwnodeIDs;

  // Multi forests: Vector of split values for each node;
  // for terminal nodes the prediction value is saved here
  std::vector<std::vector<double>> split_muwvalues;

  // Multi forests: Indicator that tells whether a split is multi-way or binary:
  std::vector<size_t> muw_inds;
  
  // Multi forests: Vector of child nodes assigned to each class:
  std::vector<std::vector<size_t>> assigned_classes;
  
  // Multi forests: Vector of classes available at each node:
  std::vector<std::vector<size_t>> classes_at_nodes;
  
  // Classes of the dependent variable and classIDs for responses
  const std::vector<double>* class_values;
  const std::vector<uint>* response_classIDs;
  const std::vector<std::vector<size_t>>* sampleIDs_per_class;

  // Class counts in terminal nodes. Empty for non-terminal nodes.
  std::vector<std::vector<double>> terminal_class_counts;

  // Splitting weights
  const std::vector<double>* class_weights;

  std::vector<size_t> counter;
  std::vector<size_t> counter_per_class;
  
  // Variable importance for multiway and binary splits:
  std::vector<double>* var_imp_multiway;
  std::vector<double>* var_imp_binary;  
  
};

} // namespace diversityForestD

#endif /* TREEPROBABILITY_H_ */
