import dendropy
from dendropy.calculate import treecompare
import math

def compare_tree_nodes(tree1, tree2, length_threshold=1e-14):
    """
    Recursively compare nodes of two trees, finding symmetric differences,
    intersections, and branch length discrepancies within a threshold.

    Parameters:
        tree1 (dendropy.Tree): First tree.
        tree2 (dendropy.Tree): Second tree.
        length_threshold (float): Threshold for acceptable branch length differences.

    Returns:
        None
    """
    def visit(node1, node2, path="root"):
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

        # Recurse into common children
        for label in common_labels:
            child1 = next((child for child in node1.child_nodes() if child.taxon and child.taxon.label == label), None)
            child2 = next((child for child in node2.child_nodes() if child.taxon and child.taxon.label == label), None)
            visit(child1, child2, path=f"{path}/{label}")

    visit(tree1.seed_node, tree2.seed_node)

def calculate_rf_distance(tree_file1, tree_file2, length_threshold=1e-14):
    """
    Calculate the Robinson-Foulds (RF) distance between two phylogenetic trees,
    (the original and the new) and compare their nodes and branch lengths.

    Parameters:
        tree_file1 (str): Path to the first .tree file (the original).
        tree_file2 (str): Path to the second .tree file (ours).
        length_threshold (float): Threshold for acceptable branch length differences.

    Returns:
        int: The RF distance between the two trees.
    """
    # Shared taxon namespace
    taxa = dendropy.TaxonNamespace()

    # Trees --> shared taxon namespace
    tree1 = dendropy.Tree.get(path=tree_file1, schema="newick", taxon_namespace=taxa)
    tree2 = dendropy.Tree.get(path=tree_file2, schema="newick", taxon_namespace=taxa)

    # Unweighted Robinson-Foulds distance
    rf_distance = treecompare.symmetric_difference(tree1, tree2)

    # Compare nodes and branch lengths for diffs
    print("Comparing node differences and branch lengths:")
    compare_tree_nodes(tree1, tree2, length_threshold=length_threshold)

    return rf_distance

# File paths of the two trees
tree_file1 = "og_output_tree.tree"
tree_file2 = "ours_output_tree.tree"

# Calculating
rf_distance = calculate_rf_distance(tree_file1, tree_file2, length_threshold=1e-14)

# RF distance (quant) output 
print(f"The Robinson-Foulds distance is: {rf_distance}")

if rf_distance == 0:
    print("The trees are identical.")
else:
    print(f"The Robinson-Foulds distance is: {rf_distance}")


