test_that("SubgridStats calculates statistics correctly", {
  TestMat <- matrix(c(1, 0, 1, 0, 1,
                       0, 1, 0, 1, 0,
                       1, 0, 1, 0, 1,
                       0, 1, 0, 1, 0,
                       1, 0, 1, 0, 1),
                     nrow = 5)

  TestRast <- terra::rast(TestMat)

  DF1 = SubgridStats(TestRast, factv1 = 5, facth1 = 5, statistic = "mean")
  DF2 = SubgridStats(TestRast, factv1 = 5, facth1 = 5, statistic = "var")
  DF3 = SubgridStats(TestRast, factv1 = 5, facth1 = 5, statistic = "sum")
  expect_equal(DF1$Mean[1], 0.52)
  expect_equal(DF2$Var[1], 0.26)
  expect_equal(DF3$Sum[1], 13)
})
