#' @export
quatrads <- function() {
  
  nodes <- CausalMST::quatrad_nodes

  nodes <- dplyr::mutate(
    dplyr::filter(
      tidyr::gather(
        dplyr::select(nodes, -cmst, -weight),
        quatrad, node, -id, -model),
      !is.na(node)),
    id = as.integer(stringr::str_replace(id, "G", "")),
    quatrad = factor(quatrad, levels = unique(quatrad)))
  
  node_id <- dplyr::select(
    tidyr::spread(
      dplyr::select(nodes, 
                    node, quatrad, id), 
      quatrad, id),
    -node)
  list(nodes = nodes, node_id = node_id)
}