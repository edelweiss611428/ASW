dx_faithful = function() stats::dist(scale(datasets::faithful))

# A small continuous data set: no duplicated points, so no exact ties in the
# ASW and the efficient and original variants must agree exactly.
dx_random = function(n = 60, seed = 1L) {
  set.seed(seed)
  stats::dist(matrix(stats::rnorm(n * 2L), ncol = 2L))
}

test_that("effOSil reproduces OSil exactly on continuous data", {
  dx = dx_random()
  for (k in 2:6) {
    fast = effOSil(dx, K = k, variant = "efficient")
    slow = effOSil(dx, K = k, variant = "original")
    expect_identical(fast$best_clustering, slow$best_clustering)
    expect_equal(fast$best_asw, slow$best_asw)
    expect_identical(fast$nIter, slow$nIter)
  }
})

test_that("scalOSil reproduces FOSil exactly", {
  dx = dx_random(80L)
  n = 24L
  for (k in 2:4) {
    set.seed(42)
    fast = scalOSil(dx, K = k, n = n, ns = 3, variant = "scalable")
    set.seed(42)
    slow = scalOSil(dx, K = k, n = n, ns = 3, variant = "original")
    expect_identical(fast$best_clustering, slow$best_clustering)
    expect_equal(fast$best_asw, slow$best_asw)
  }
})

test_that("asw() agrees with cluster::silhouette()", {
  skip_if_not_installed("cluster")
  dx = dx_faithful()
  for (k in 2:5) {
    cl = effOSil(dx, K = k)$best_clustering
    expect_equal(asw(cl, dx), mean(cluster::silhouette(cl, dx)[, 3]))
  }
})

test_that("Silhouette() matches cluster::silhouette()", {
  skip_if_not_installed("cluster")
  dx = dx_faithful()
  cl = effOSil(dx, K = 3)$best_clustering
  sw = Silhouette(cl, dx)
  ref = cluster::silhouette(cl, dx)

  expect_s3_class(sw, "silhouette")
  expect_equal(mean(sw[, "sil_width"]), asw(cl, dx))
  expect_equal(as.numeric(sw[, "sil_width"]), as.numeric(ref[, "sil_width"]))
  expect_equal(as.numeric(sw[, "neighbor"]), as.numeric(ref[, "neighbor"]))
  expect_true(all(sw[, "neighbor"] != sw[, "cluster"]))
})

test_that("the reported ASW matches the reported clustering", {
  dx = dx_random(50L)
  for (fit in list(effOSil(dx, K = 2:5), PAMSil(dx, K = 2:5))) {
    expect_equal(fit$best_asw, asw(fit$best_clustering, dx))
    for (i in seq_along(fit$asw)) {
      expect_equal(unname(fit$asw[i]), asw(fit$clusterings[, i], dx))
    }
  }
})

test_that("scalOSil is consistent when rep > 1", {
  dx = dx_random(70L)
  set.seed(7)
  fit = scalOSil(dx, K = 3, n = 20L, ns = 3, rep = 4)
  expect_equal(fit$best_asw, asw(fit$best_clustering, dx))
})

test_that("degenerate inputs are handled", {
  # every point coincident: all silhouette widths are 0/0, defined as 0
  dx = stats::dist(matrix(0, nrow = 8L, ncol = 2L))
  expect_equal(asw(rep(1:2, each = 4L), dx), 0)
  expect_true(is.finite(effOSil(dx, K = 2)$best_asw))

  # singleton clusters get width 0, but their neighbour is still reported
  dxr = dx_random(6L)
  sw = Silhouette(c(1L, 2L, 2L, 2L, 3L, 2L), dxr)
  expect_equal(unname(sw[1L, "sil_width"]), 0)
  expect_true(sw[1L, "neighbor"] %in% c(2L, 3L))
})

test_that("invalid arguments are rejected", {
  dx = dx_random(30L)

  expect_error(effOSil(as.matrix(dx)), "dist")
  expect_error(effOSil(dx, K = 1), "at least 2")
  expect_error(effOSil(dx, K = 1000), "cannot exceed")
  expect_error(effOSil(dx, K = c(2, 2)), "duplicated")
  expect_error(effOSil(dx, K = numeric(0)), "non-empty")
  expect_error(effOSil(dx, variant = "nope"))

  expect_error(scalOSil(dx, n = 1), "at least 2")
  expect_error(scalOSil(dx, n = 1000), "cannot exceed")
  expect_error(scalOSil(dx, K = 20, n = 10), "subsample size")
  expect_error(scalOSil(dx, ns = 0), "at least 1")

  bad = dx
  bad[1L] = NA_real_
  expect_error(effOSil(bad), "missing")

  bad2 = dx
  bad2[1L] = -1
  expect_error(effOSil(bad2), "negative")
})

test_that("asw() accepts arbitrary cluster labels", {
  dx = dx_random(40L)
  cl = effOSil(dx, K = 3)$best_clustering
  expect_equal(asw(letters[cl], dx), asw(cl, dx))
})

test_that("Init picks the method with the highest ASW", {
  dx = dx_random(50L)
  methods = c("pam", "average", "complete", "single")
  fit = Init(dx, 3, methods)

  expect_named(fit$all_asw, methods)
  expect_equal(fit$asw, max(fit$all_asw))
  expect_identical(fit$method, methods[which.max(fit$all_asw)])
  expect_equal(fit$asw, asw(fit$clustering, dx))
  expect_length(unique(fit$clustering), 3L)

  expect_error(Init(dx, 3, "nope"), "Unsupported")
  expect_error(Init(dx, 1), "at least 2")
  expect_error(Init(dx, 1000), "cannot exceed")
})

test_that("print returns its argument invisibly", {
  fit = effOSil(dx_random(30L), K = 2:3)
  expect_output(print(fit), "effOSil")
  expect_identical(withVisible(print(fit))$visible, FALSE)
})
