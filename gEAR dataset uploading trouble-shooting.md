gEAR dataset uploading trouble-shooting


1. "can't complete the uploading"
* Make sure there are no rows without ensemble ID in either expression or genes tab/sheet
* Do not include the greek letters in the metadata sheet.
 


 2. Can’t curate umap_static by group
* Make sure there are no special characters in the column, for example “/”, “+”, space.
This dataset failed due to “ IP/IB,Oc90+,Immune cells”
* Make sure the cluster names are not in numeric format


  



3. In order to enable the “primary analysis” function of single cell workbench, one of the column names MUST be in this list. ['cluster', 'cell_type', 'cluster_label', 'subclass_label', 'joint_cluster_round4_annot']. As long as one of the tag is found in this list, the primary analysis will use that tag to display the cluster to compare with.


4. Don’t put “,” or “ “ in the column name, otherwise, the multi-genes display function will have problems.