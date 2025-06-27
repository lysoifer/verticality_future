reproj_df = function(df, from, to, method, torast) {
  d = df %>% 
    rast(crs = from) %>%
    project(to, method = method)
  
  if(!torast) {d = as.data.frame(d, xy = T)}
  
  return(d)
}
