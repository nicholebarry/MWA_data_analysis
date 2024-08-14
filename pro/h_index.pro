pro h_index

  ;Last updated 2-9-2019 from ADS citation counts

  paper_list = ['barry_2016','barry_2019','trott_2019b','wilensky_2019','byrne_2019','morales_2019','trott_2019a',$
    'kerrigan_2018','Li_2018','lenc_2017','kapinska_2017','sourabh_2016','beardsley_2016','carroll_2016',$
    'lenc_2016','ewall-wice_2016','jacobs_2016','offringa_2016','pober_2016','trott_2016','thyagarajan_2015',$
    'dillon_2015','thyagarajan_2015b','offringa_2015','seidler_2014','mock_2014','szydagis_2011'];'barry_thesis'

  citation_list = [56,4,1,1,8,5,3,7,10,22,10,5,72,15,44,38,49,39,47,61,49,75,79,73,21,21,61];3

  print, strtrim(N_elements(paper_list),2) + ' papers with ' + strtrim(total(citation_list),2) + ' total citations'
  print, [transpose(paper_list),string(transpose(citation_list))]

  sorted_citation_list = sort(citation_list)
  index_num = indgen(N_elements(citation_list))+1

  h_index_num=0
  ind = 0
  while h_index_num EQ 0 do begin
    if (citation_list[reverse(sorted_citation_list)])[ind] LE index_num[ind] then h_index_num = index_num[ind]
    ind=ind+1
  endwhile

print, 'h-index is ' + strtrim(h_index_num,2)

end