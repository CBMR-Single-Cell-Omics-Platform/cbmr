#' Download metadata from labguru
#'
#' Download metadata from labguru when given a experiment id or
#' an experiment name
#'
#' @param expId numeric experimental id, fastest way to get data
#' @param expName character of experimental name
#' @export
#' @examples
#' \dontrun{
#' labguru_set_token(token = "abcdefg",
#' server = "https://sund.labguru.com")
#' getMetadata(expId = 1)
#' }

getMetadata <- function(expId = NULL, expName = NULL){
  if (is.null(expId)){
    if (is.null(expName)) stop("expId or expName must be set.")
    allExps <- LabguruR::labguru_list_experiments()
    expId <- allExps[allExps[["title"]] == expName, "id"]
    if (length(expId) == 0) {
      stop("expName matched no experiments. Check name is correct, or supply expId instead")
    }
    if (length(expId) > 1) {
      stop("expName matched multiple experiments. Supply expId instead")
    }
  }

  procedures <- LabguruR::labguru_get_experiment(expId)$experiment_procedures
  procedureId <- procedures[
    procedures[["experiment_procedure.name"]] == "Procedure",
    "experiment_procedure.id",
    drop = TRUE]

  # this function is not exported from the labguru package, but the package has
  # not been updated to match the current API. When the package updates
  # labguru_get_experiment_procedure() will replace this function
  labguru_get_by_id <- function (type,
                                 id,
                                 server = Sys.getenv("LABGURU_SERVER"),
                                 token = Sys.getenv("LABGURU_TOKEN"))
  {
    base_url <- server
    path <- paste0("/api/v1/", type, "/", id)
    query <- paste0("token=", token)
    url <- httr::modify_url(url = base_url, path = path, query = query)
    url <- utils::URLencode(url)
    resp <- httr::GET(url)
    if (httr::http_type(resp) != "application/json") {
      stop("API did not return JSON", call. = FALSE)
    }
    parsed <- jsonlite::fromJSON(httr::content(resp, as = "text"),
                                 simplifyVector = FALSE,
                                 simplifyDataFrame = TRUE,
                                 flatten = TRUE)
    if (httr::http_error(resp)) {
      stop(sprintf("API request failed [%s]\n%s", parsed$status,
                   parsed$error), call. = FALSE)
    }
    parsed
  }

  parsed <- labguru_get_by_id(type = "sections",
                              id = procedureId)$elements
  sheetId <- parsed[parsed[["element_type"]] == "excel", "id", drop = TRUE]

  datasheetJSON <- LabguruR::labguru_get_element(sheetId)$data

  datasheetList <- jsonlite::fromJSON(
    jsonlite::fromJSON(txt = datasheetJSON)$spread
  )
  datasheetList <- datasheetList$sheets$Sheet1$data$dataTable
  getValues <- function(x){
    lapply(x, magrittr::extract2, "value")
  }
  filterFn <- function(x) !is.null(x$`1`)

  datasheet <- lapply(datasheetList, getValues)
  datasheet <- Filter(f = filterFn, datasheet)

  out <- as.data.frame(matrix(nrow = length(datasheet) - 1,
                              ncol = length(datasheet$`0`)))
  colnames(out) <- unlist(datasheet$`0`)
  datasheet <- datasheet[-1]

  for (i in seq_along(datasheet)){
    for (j in names(datasheet[[i]])){
      idx <- as.integer(j) + 1
      val <- datasheet[[i]][[j]]
      if(is.null(val)) val <- NA
      out[i, idx] <- val
    }
  }
  data.table::setDT(out)
  return(out)
}
