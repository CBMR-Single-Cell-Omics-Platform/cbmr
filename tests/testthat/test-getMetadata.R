test_that("download example metadata", {
  expectedResult <- structure(list(
    `Sample ID` = c("Ctr1", "Ctrl2", "Ctrl3", "Ctr4", "Ctrl5", "WT1", "WT2",
                    "WT3", "WT4", "WT5", "WT6", "WT7"),
    `Basespace ID` = c("101_1", "101_2", "101_3", "101_4", "101_5", "101_6",
                       "101_7", "101_8", "101_9", "101_10", "101_11", "101_12"),
    Species = c("Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse",
                "Mouse", "Mouse", "Mouse", "Mouse", "Mouse"),
    Project = c("NADcap", "NADcap", "NADcap", "NADcap",
                "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT"),
    `Nucleic acid` = c("DNA", "DNA", "DNA", "DNA", "DNA", "DNA", "DNA", "DNA",
                       "DNA", "DNA", "DNA", "DNA"),
    `Extraction protocol` = c("RNeasy", "RNeasy", "RNeasy", "RNeasy", "RNeasy",
                              "RNeasy", "RNeasy", "RNeasy", "RNeasy", "RNeasy",
                              "RNeasy", "RNeasy"),
    `Condition 1` = c("Ctrl", "Ctrl", "Ctrl", "Ctrl", "Ctrl", "WT", "WT", "WT",
                      "WT", "WT", "WT", "WT"),
    `Condition 2` = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)),
    row.names = c(NA, -12L), class = c("data.frame"))
  data.table::setDT(expectedResult)

  x <- getMetadata(35)

  expect_equal(x, expectedResult)
})
