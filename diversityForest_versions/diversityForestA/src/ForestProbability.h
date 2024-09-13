/*-------------------------------------------------------------------------------
 This file is part of diversityForestA.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of divfor is distributed under MIT license and the
 R package "diversityForestA" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#ifndef FORESTPROBABILITY_H_
#define FORESTPROBABILITY_H_

#include <map>
#include <utility>
#include <vector>

#include "globals.h"
#include "Forest.h"
#include "TreeProbability.h"

namespace diversityForestA {

class ForestProbability: public Forest {
public:
  ForestProbability() = default;

  ForestProbability(const ForestProbability&) = delete;
  ForestProbability& operator=(const ForestProbability&) = delete;

  virtual ~ForestProbability() override = default;

  void loadForest(size_t dependent_varID, size_t num_trees,
      std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
      std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
	  std::vector<std::vector<size_t>>& forest_split_types, std::vector<std::vector<std::vector<size_t>>>& forest_split_multvarIDs, 
	  std::vector<std::vector<std::vector<std::vector<bool>>>>& forest_split_directs, 
	  std::vector<std::vector<std::vector<std::vector<double>>>>& forest_split_multvalues, 
	  std::vector<std::vector<std::vector<size_t>>>& forest_child_muwnodeIDs,
	  std::vector<std::vector<std::vector<double>>>& forest_split_muwvalues, 	
      std::vector<double>& class_values, std::vector<std::vector<std::vector<double>>>& forest_terminal_class_counts,
      std::vector<bool>& is_ordered_variable);

  std::vector<std::vector<std::vector<double>>> getTerminalClassCounts() const;

  const std::vector<double>& getClassValues() const {
    return class_values;
  }

  void setClassWeights(std::vector<double>& class_weights) {
    this->class_weights = class_weights;
  }
  
  std::vector<std::vector<std::vector<double>>> getSplitMuwValues() {
    std::vector<std::vector<std::vector<double>>> result;
    for (auto &tree : trees) {
      TreeProbability *treeProb = dynamic_cast<TreeProbability *>(tree.get());
      if (treeProb)
      {
        result.push_back(treeProb->getSplitMuwValues());
      }
    }
    return result;
  }

  std::vector<std::vector<std::vector<size_t>>> getChildMuwNodeIDs() {
    std::vector<std::vector<std::vector<size_t>>> result;
    for (auto &tree : trees) {
      TreeProbability *treeProb = dynamic_cast<TreeProbability *>(tree.get());
      if (treeProb)
      {
        result.push_back(treeProb->getChildMuwNodeIDs());
      }
    }
    return result;
  }
  
  std::vector<std::vector<size_t>> getMuwInds() {
    std::vector<std::vector<size_t>> result;
    for (auto &tree : trees) {
      TreeProbability *treeProb = dynamic_cast<TreeProbability *>(tree.get());
      if (treeProb)
      {
        result.push_back(treeProb->getMuwInds());
      }
    }
    return result;
  }

  const std::vector<double>& getVariableImportanceMuwMultiway() const {
    return var_imp_multiway;
  }
  
  const std::vector<double>& getVariableImportanceMuwBinary() const {
    return var_imp_binary;
  }

protected:
  void initInternal(std::string status_variable_name) override;
  void growInternal() override;
  void allocatePredictMemory() override;
  void predictInternal(size_t sample_idx) override;
  void computePredictionErrorInternal() override;
  void saveToFileInternal(std::ofstream& outfile) override;

  void computeImportanceMuw() override;

  void computeTreeImportanceMuwInThread(uint thread_idx, std::vector<double>& importance_multiway,
     std::vector<double>& importance_binary);

  // Classes of the dependent variable and classIDs for responses
  std::vector<double> class_values;
  std::vector<uint> response_classIDs;
  std::vector<std::vector<size_t>> sampleIDs_per_class;

  // Splitting weights
  std::vector<double> class_weights;

  // Variable importance for multiway and binary splits:
  std::vector<double> var_imp_multiway;
  std::vector<double> var_imp_binary;

private:
  const std::vector<double>& getTreePrediction(size_t tree_idx, size_t sample_idx) const;
  size_t getTreePredictionTerminalNodeID(size_t tree_idx, size_t sample_idx) const;
};

} // namespace diversityForestA

#endif /* FORESTPROBABILITY_H_ */
