test_that("Moran's I works", {
  TestMat <- matrix(c(1, 0, 1, 0, 1,
                      0, 1, 0, 1, 0,
                      1, 0, 1, 0, 1,
                      0, 1, 0, 1, 0,
                      1, 0, 1, 0, 1),
                    nrow =5)
  MoransI(TestMat, r1 = 1)
  expect_equal(MoransI(TestMat, r1 = 1), -1)
})
