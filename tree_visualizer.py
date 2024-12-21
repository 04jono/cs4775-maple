# from ete3 import Tree, TreeStyle

# # Load the tree file
# tree = Tree("ogMAPLE_output_tree.tree", format=1)

# # Configure a circular style
# ts = TreeStyle()
# ts.mode = "c"  # Circular mode
# ts.show_leaf_name = True  # Show leaf names
# ts.scale = 10  # Adjust scale for large trees

# # Render the tree
# tree.show(tree_style=ts)

'''picking 20 random individuals'''
# import random
# from ete3 import Tree, TreeStyle, NodeStyle
# import re
# from matplotlib import colormaps
# import numpy as np

# # Step 1: Load the tree
# tree_file = "ogMAPLE_output_tree.tree"
# ete_tree = Tree(tree_file, format=1)

# # Step 2: Randomly sample 20 leaf nodes
# all_leaves = [leaf for leaf in ete_tree]
# sampled_leaves = random.sample(all_leaves, 20)

# # Save the names of the sampled individuals
# sampled_individuals = [leaf.name for leaf in sampled_leaves]
# print("Sampled Individuals:")
# print(sampled_individuals)

# # Step 3: Prune the tree to include only the sampled individuals
# pruned_tree = ete_tree.copy()
# pruned_tree.prune(sampled_individuals, preserve_branch_length=True)

# # Step 4: Map individuals to distinct colors for clarity
# colors = colormaps["tab20"]
# num_individuals = len(sampled_individuals)
# individual_colors = {
#     sampled_individuals[i]: tuple((np.array(colors(i / num_individuals)[:3]) * 255).astype(int))
#     for i in range(num_individuals)
# }

# # Step 5: Function to extract the individual code (EPI_ISL_)
# def extract_individual_code(name):
#     match = re.search(r"EPI_ISL_\d+", name)
#     return match.group(0) if match else name  # Return the EPI_ISL code or the full name if not found

# # Step 6: Annotate the tree with colors and remove dots at nodes, show short names
# def set_node_style(node, color):
#     style = NodeStyle()
#     style["fgcolor"] = f"rgb({color[0]},{color[1]},{color[2]})"  # Convert to RGB string
#     style["size"] = 0  # Remove node dot by setting size to 0
#     node.set_style(style)
    
#     # Extract the individual code for the label
#     short_name = extract_individual_code(node.name)
#     node.name = short_name  # Replace the full name with the short code

# for leaf in pruned_tree:
#     if leaf.name in individual_colors:
#         set_node_style(leaf, individual_colors[leaf.name])

# # Step 7: Visualize and save the tree
# ts = TreeStyle()
# ts.mode = "c"  # Circular mode
# ts.show_leaf_name = True
# ts.scale = 50  # Adjust scale for visibility

# # Save the tree as an image
# output_file = "pruned_ogphylogenetic_tree.png"  # Specify the output file name
# pruned_tree.render(output_file, tree_style=ts)

# print(f"Tree saved as {output_file}")

'''with the 20 individuals picked here is the final'''
import re
from ete3 import Tree, TreeStyle, NodeStyle
from matplotlib import colormaps
import numpy as np

# Step 1: Load the tree
tree_file = "ogMAPLE_output_tree.tree"
ete_tree = Tree(tree_file, format=1)

# Step 2: Specify the list of sampled individuals directly
sampled_individuals = [
    'hCoV-19/China/FU-P242-2/2020|EPI_ISL_19023709|2020-02-03', 
    'hCoV-19/Shanghai/SH0026/2020|EPI_ISL_416335|2020-02-02', 
    'hCoV-19/Wuhan/0125-A153/2020|EPI_ISL_493154|2020-01-25', 
    'hCoV-19/Wuhan/0125-A167/2020|EPI_ISL_493158|2020-01-25', 
    'hCoV-19/Zhejiang/YW03_Q1/2020|EPI_ISL_2433283|2020-02-02', 
    'hCoV-19/Sichuan/SC-YB-085/2020|EPI_ISL_451398|2020-01-31', 
    'hCoV-19/Wuhan/0125-A148/2020|EPI_ISL_493152|2020-01-25', 
    'hCoV-19/Zhejiang/HZ103/2020|EPI_ISL_422425|2020-01-24', 
    'hCoV-19/Wuhan/IPBCAMS-WH-01/2019|EPI_ISL_402123|2019-12-24', 
    'hCoV-19/Sichuan/SC-WCH3-259/2020|EPI_ISL_451385|2020-02-03', 
    'hCoV-19/Shandong/LY006-2/2020|EPI_ISL_414939|2020-01-25', 
    'hCoV-19/Anhui/226/2020|EPI_ISL_1069216|2020-02-04', 
    'hCoV-19/Hangzhou/HZCDC0091L/2020|EPI_ISL_421227|2020-01-21', 
    'hCoV-19/Guangzhou/GZMU0036/2020|EPI_ISL_429104|2020-02-01', 
    'hCoV-19/Hangzhou/ZJU-08/2020|EPI_ISL_416473|2020-01-26', 
    'hCoV-19/Guangdong/20SF014/2020|EPI_ISL_403934|2020-01-15', 
    'hCoV-19/Shandong/2020C1240116/2020|EPI_ISL_962875|2020-02-02', 
    'hCoV-19/Wuhan/0126-C7/2020|EPI_ISL_493168|2020-01-26', 
    'hCoV-19/Guangdong/SZTH-002/2020|EPI_ISL_406593|2020-01-13', 
    'hCoV-19/Hangzhou/ZJU-01/2020|EPI_ISL_415709|2020-01-25'
]

# Step 3: Prune the tree to include only the sampled individuals
pruned_tree = ete_tree.copy()
pruned_tree.prune(sampled_individuals, preserve_branch_length=True)

# Step 4: Map individuals to distinct colors for clarity
colors = colormaps["tab20"]
num_individuals = len(sampled_individuals)
individual_colors = {
    sampled_individuals[i]: tuple((np.array(colors(i / num_individuals)[:3]) * 255).astype(int))
    for i in range(num_individuals)
}

# Step 5: Function to extract the individual code (EPI_ISL_)
def extract_individual_code(name):
    match = re.search(r"EPI_ISL_\d+", name)
    return match.group(0) if match else name  # Return the EPI_ISL code or the full name if not found

# Step 6: Annotate the tree with colors and remove dots at nodes, show short names
def set_node_style(node, color):
    style = NodeStyle()
    style["fgcolor"] = f"rgb({color[0]},{color[1]},{color[2]})"  # Convert to RGB string
    style["size"] = 0  # Remove node dot by setting size to 0
    node.set_style(style)
    
    # Extract the individual code for the label
    short_name = extract_individual_code(node.name)
    node.name = short_name  # Replace the full name with the short code

for leaf in pruned_tree:
    if leaf.name in individual_colors:
        set_node_style(leaf, individual_colors[leaf.name])

# Step 7: Visualize and save the tree
ts = TreeStyle()
ts.mode = "c"  # Circular mode
ts.show_leaf_name = True
ts.scale = 50  # Adjust scale for visibility

# Save the tree as an image
output_file = "pruned_og_final.png"  # Specify the output file name
pruned_tree.render(output_file, tree_style=ts)

print(f"Tree saved as {output_file}")

