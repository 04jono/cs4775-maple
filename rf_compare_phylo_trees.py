import dendropy
from dendropy.calculate import treecompare
import math

def compare_tree_nodes(tree1, tree2, length_threshold=1e-14):
    """
    Recursively compare nodes of two trees, finding symmetric differences,
    intersections, and branch length discrepancies within a threshold.
    """
    total_length_diff = 0  

    def visit(node1, node2, path="root"):
        nonlocal total_length_diff
        if node1 is None or node2 is None:
            print(f"Mismatch at {path}: One node is missing.")
            return

        # Child labels and branch lengths
        children1 = {child.taxon.label for child in node1.child_nodes() if child.taxon}
        children2 = {child.taxon.label for child in node2.child_nodes() if child.taxon}
        
        lengths1 = {child.taxon.label: child.edge_length for child in node1.child_nodes() if child.taxon}
        lengths2 = {child.taxon.label: child.edge_length for child in node2.child_nodes() if child.taxon}

        # Comparing child labels
        print(f"Node: {path}")
        print(f"  Only in tree1: {children1 - children2}")
        print(f"  Only in tree2: {children2 - children1}")
        print(f"  Common: {children1.intersection(children2)}")

        # Branch lengths common children
        common_labels = children1.intersection(children2)
        for label in common_labels:
            length1 = lengths1[label]
            length2 = lengths2[label]
            if length1 is not None and length2 is not None:
                if not math.isclose(length1, length2, rel_tol=0, abs_tol=length_threshold):
                    print(f"  Branch length discrepancy for {label}: Tree1={length1}, Tree2={length2}")
                    total_length_diff += abs(length1 - length2)

        # Recurse into common children
        for label in common_labels:
            child1 = next((child for child in node1.child_nodes() if child.taxon and child.taxon.label == label), None)
            child2 = next((child for child in node2.child_nodes() if child.taxon and child.taxon.label == label), None)
            visit(child1, child2, path=f"{path}/{label}")

    visit(tree1.seed_node, tree2.seed_node)

    return total_length_diff

def calculate_rf_and_rfl_distance(tree_file1, tree_file2, length_threshold=1e-14):
    """
    Calculate Robinson-Foulds (RF), maximum RF, and Branch Length Robinson-Foulds (RFL) distances
    between two phylogenetic trees, along with node and branch length comparisons.
    """
    # Shared taxon namespace
    taxa = dendropy.TaxonNamespace()

    # Trees --> shared taxon namespace
    tree1 = dendropy.Tree.get(path=tree_file1, schema="newick", taxon_namespace=taxa)
    tree2 = dendropy.Tree.get(path=tree_file2, schema="newick", taxon_namespace=taxa)

    # Unweighted Robinson-Foulds distance (RF)
    rf_distance = treecompare.symmetric_difference(tree1, tree2)

    # Maximum RF distance (2(n - 3))
    n_taxa = len(taxa)
    max_rf_distance = 2 * (n_taxa - 3)

    print("Comparing node differences and branch lengths:")
    total_length_diff = compare_tree_nodes(tree1, tree2, length_threshold=length_threshold)

    rfl_distance = rf_distance + total_length_diff

    return rf_distance, max_rf_distance, rfl_distance

# Input tree files - can be replaced with any two .tree files
tree_file1 = "_tree.tree"  # Our tree
tree_file2 = "may12021.fasta.treefile"  # IQ-TREE

rf_distance, max_rf_distance, rfl_distance = calculate_rf_and_rfl_distance(tree_file1, tree_file2, length_threshold=1e-14)

print(f"The Robinson-Foulds distance (RF) is: {rf_distance}")
print(f"The Maximum Robinson-Foulds distance (Max RF) is: {max_rf_distance}")
print(f"The Branch Length Robinson-Foulds distance (RFL) is: {rfl_distance}")

if rf_distance == 0 and rfl_distance == 0:
    print("The trees are identical.")
else:
    print(f"RF Distance: {rf_distance}, Max RF Distance: {max_rf_distance}, RFL Distance: {rfl_distance}")
