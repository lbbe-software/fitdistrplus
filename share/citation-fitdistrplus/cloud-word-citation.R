
library(RefManageR)

mescitations <- RefManageR::ReadBib("citations_fitdistrplus_2025-2026-1.bib")

mescitations <- as.data.frame(mescitations)
mescitations$title <- tolower(mescitations$title)
mescitations$journal <- tolower(mescitations$journal)

head(sort(table(mescitations$journal), decreasing =TRUE), 10)

require(wordcloud)
require(tm)
makemycloud <- function(data, colname,  min.freq, ...)
{
  data <- data[, colname]
  #print(head(data))
  docs <- Corpus(VectorSource(data)) 
  #print(head(docs))
  df <- getfreq(docs)
  forbidden <- c("for", "and", "the", "a", "with", "from", "using", "via", "des",
                 "The", "based", "'de", "also", "pour", "les", "new", "across",
                 "between", "during", "through", "under", "new", "how", "but", "not",
                 "are", "via", "after", "into", "its", "within", "journal")
  df <- subset(df,! word %in% forbidden)
  
  print(head(df))
  
  n <- sum(df$freq >= min.freq)
  wordcloud(words=df$word, freq=df$freq, colors=rainbow(n), random.color=TRUE,
            min.freq=min.freq, max.words=Inf, ...) 
  
}

getfreq <- function(doc)
{
  #minuscule
  #docs <- tm_map(docs, content_transformer(tolower)) 
  #term matrix
  dtm <- TermDocumentMatrix(doc, control=list(tolower=FALSE))
  dtm <- as.matrix(dtm)
  
  words <- sort(rowSums(dtm),decreasing=TRUE)
  df <- data.frame(word = names(words),freq=words) 
  df
}

getfreq(Corpus(VectorSource(mescitations$title)) )

png("citations_fitdistrplus_2025-2026_title.png", "500", "500")
makemycloud(mescitations, "title", min.freq=4, fixed.asp=FALSE, rot.per=0)
dev.off()


png("citations_fitdistrplus_2025-2026_journal.png", "500", "500")
makemycloud(mescitations, "journal", min.freq=2, fixed.asp=FALSE, rot.per=0)
dev.off()
