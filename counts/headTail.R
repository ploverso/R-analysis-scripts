# Some simple functions for more succiently returning heads and tails of strings

strtail <- function(s, n) {
	if(n < 0){
		return(substring(s, 1 - n))
	}
	else{
		return(substring(s, nchar(s) - n + 1))
	}
}
strhead <- function(s, n) {
	if(n < 0){
		return(substr(s, 1, nchar(s) + n))
	}
	else{
		return(substr(s, 1, n))
	}
}
