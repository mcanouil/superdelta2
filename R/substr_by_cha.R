#' Divide character strings given a certain character.
#'
#' This is a handy function to divide character strings given a certain character. This is especially useful
#' when sometimes gene names show up in the form of "GENESYMBOL|ENTREZID"
#'
#' @param cha_str A character string or a vector of character string.
#' @param delim Delimiter character. Default is "|".
#' @param type Return the partial character string before or after the delimiter.
#'
#' @return Return the desired substring (an individual or a vector).
#'
#' @export
#'
#' @keywords utilities misc
#'
#' @examples
#' substr_by_cha("ABCB1|5243", delim = "|", type = "before")
#' substr_by_cha(c("ABCB1-5243", "ABL1-25", "ABL2-27"), delim = "-", type = "after")
substr_by_cha <- function(cha_str, delim = "|", type = c("before", "after")) {
  if (!is.character(cha_str)) {
    stop("Input must be a character vector!")
  }
  cha_to_return <- character()
  if (length(cha_str) != 0) {
    # cha_to_return <- character()
    for (i in seq_along(cha_str)) {
      loc <- which(strsplit(cha_str[i], "")[[1]] == delim)
      if (type == "before") {
        cha_to_return[i] <- substr(cha_str[i], 1, loc - 1)
      } else if (type == "after") {
        cha_to_return[i] <- substr(cha_str[i], loc + 1, nchar(cha_str[i]))
      } else {
        stop("Substring type needs to be corrected!")
      }
    }
  }
  return(cha_to_return)
}
