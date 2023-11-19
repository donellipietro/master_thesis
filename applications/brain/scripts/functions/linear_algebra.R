# Linear algebra ----
# |||||||||||||||||||

rodrigues_rotation_formula <- function(w, theta) {
  # Normalize the rotation axis
  w <- w / sqrt(sum(w^2))
  # Skew-symmetric matrix K
  K <- matrix(c(0, -w[3], w[2], w[3], 0, -w[1], -w[2], w[1], 0), nrow = 3)
  # Rodrigues' rotation formula for the rotation matrix
  R <- diag(3) + sin(theta) * K + (1 - cos(theta)) * K %*% K
  return(R)
}

rotation_matrix <- function(v, e = c(1, 0, 0)){

  # Normalization
  v <- v/sqrt(sum(v^2))
  
  # Calculate the cross product to find the axis of rotation
  w <- cross(e, v)
  
  # Calculate the angle of rotation (in radians)
  theta <- acos(sum(e * v))
  
  # Rotation matrix
  R <- rodrigues_rotation_formula(w, theta)
  
  return(R)
}
