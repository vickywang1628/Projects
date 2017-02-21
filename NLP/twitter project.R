rm(list=ls())
library("twitteR")
library("ROAuth")
library('tm')
library('lda')


download.file(url="http://curl.haxx.se/ca/cacert.pem", destfile="cacert.pem")
# reqURL <- "https://api.twitter.com/oauth/request_token"
# accessURL <- "http://api.twitter.com/oauth/access_token"
# authURL <- "http://api.twitter.com/oauth/authorize"
# consumerKey <- "gL6E6Jc6DRwEcnvaGVxMXw"
# consumerSecret <- "H98xsByVPXMUlTaKnOUNll0YXjjmY7d2g3TkaAg2dE"
# twitCred <- OAuthFactory$new(consumerKey=consumerKey,
#                              consumerSecret=consumerSecret,
#                              requestURL=reqURL,
#                              accessURL=accessURL,
#                              authURL=authURL)
# twitCred$handshake(cainfo="cacert.pem")
# registerTwitterOAuth(twitCred)

load(file="twitter authentication.Rdata")
registerTwitterOAuth(twitCred)
setup_twitter_oauth(twitCred)



ut <- userTimeline('barackobama', n=1000,cainfo="cacert.pem")
r_stats<- searchTwitter("#Rstats", n=1500, cainfo="cacert.pem")

#build a data frame
DF <- twListToDF(ut)
#Clean up the text remove punctuation and stop words
DF$text <- gsub(pattern="obama|president|http.*"
                , replacement="",x=DF$text,ignore.case=T)
#DF$text <- gsub(pattern="RT",replacement="",x=DF$text,ignore.case=T)
DF$text <- removePunctuation(DF$text)
DF$text <- removeWords(DF$text, stopwords("english"))
#build the LDA document and vocab
lex.output <- lexicalize(DF$text)
#Remove less freqeuent words
to.keep <- lex.output$vocab[word.counts(lex.output$documents, lex.output$vocab) >= 2]
#gibbs sampling model
result <- lda.collapsed.gibbs.sampler(lex.output$documents,K=10,lex.output$vocab
                                      ,num.iterations=25,0.1,0.1)

top.words <- top.topic.words(result$topics, 5, by.score=TRUE)
top.words
