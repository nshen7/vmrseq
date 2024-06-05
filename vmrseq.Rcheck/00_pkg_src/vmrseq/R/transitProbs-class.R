### =========================================================================
### transitProbs objects
### -------------------------------------------------------------------------

methods::setClass("transitProbs",
                  representation(max_dist_bp = "numeric",
                                 buffer_bp = "numeric",
                                 transit_probs = "data.frame",
                                 buffer_probs = "data.frame",
                                 train = "data.frame")
                  )

