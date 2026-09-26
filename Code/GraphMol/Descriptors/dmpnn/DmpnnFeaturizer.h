//
//  C++ port of dmpnn-cpp-fp16/tools/featurize.py: SMILES -> the graph record the
//  D-MPNN fp16 regressor consumes. The trained input format is the authority;
//  every choice here mirrors that script literally (see the notes in the .cpp).
//
#ifndef RD_DMPNN_FEATURIZER_H
#define RD_DMPNN_FEATURIZER_H

#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>

#include <cstdint>
#include <string>
#include <vector>

namespace RDKit {
namespace DmpnnFeaturizer {

constexpr unsigned int NODE_DIM = 141;       // 128-bin MMFF type one-hot ++ 13 numeric
constexpr unsigned int EDGE_DIM = 6;
constexpr unsigned int WHOLE_MOL_DIM = 11;   // on the hydrogen-explicit molecule
constexpr unsigned int DESCRIPTOR_DIM = 217; // RDKit Descriptors._descList, arcsinh'd
constexpr unsigned int MOL_DIM = WHOLE_MOL_DIM + DESCRIPTOR_DIM;

struct Graph {
  int nNodes = 0;
  int nEdges = 0;
  std::vector<float> x;         // nNodes * NODE_DIM
  std::vector<float> edgeAttr;  // nEdges * EDGE_DIM
  std::vector<std::int32_t> src, dst, rev;
  std::vector<float> mol;       // MOL_DIM
};

//! Featurize `mol` (as parsed from SMILES, no explicit hydrogens).
RDKIT_DESCRIPTORS_EXPORT Graph featurize(const ROMol &mol);

//! The graph in `predict`'s little-endian layout:
//! int32 n_nodes, n_edges, node_dim, edge_dim, mol_dim; float32 x, edge_attr;
//! int32 src, dst, rev; float32 mol.
RDKIT_DESCRIPTORS_EXPORT std::string graphBytes(const Graph &graph);

//! RDKit's 217 descriptors (Descriptors._descList order, raw values) with the
//! definition that changed after RDKit 2025.09 restored, so values match what a
//! 2025.09.x Python produced: NumHAcceptors uses the 2025 SMARTS (RDKit #9060
//! later excluded degree-3 aromatic N). NumRotatableBonds needs no change here:
//! on the hydrogen-suppressed molecule both releases agree.
RDKIT_DESCRIPTORS_EXPORT std::vector<double> rdkit217Compat2025(const ROMol &mol);

//! The 228-value molecule block alone (11 whole-molecule columns ++ 217 descriptors).
RDKIT_DESCRIPTORS_EXPORT std::vector<float> molBlock(const ROMol &mol);

}  // namespace DmpnnFeaturizer
}  // namespace RDKit

#endif
