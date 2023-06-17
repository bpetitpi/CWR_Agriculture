
### Assemble the main figure
library(magick)
a<-image_read("F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/nsp_obs_raster.png") %>% 
  image_annotate("a) ", size = 60, location = "+30+30") 
b<-image_read("F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/nsp_OR10.png") %>% 
  image_annotate("b) ", size = 60, location = "+30+30") 
c<-image_read("F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/fig_sampling_effort.png") %>% 
  image_annotate("c) ", size = 60, location = "+30+30") 
d<-image_read("F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/deficit_OR10.png") %>% 
  image_annotate("d) ", size = 60, location = "+30+30") %>% 
  image_annotate("-",size = 60, location = "+1905+245") %>% 
  image_annotate("-",size = 60, location = "+1905+410") %>%
  image_annotate("-",size = 60, location = "+1905+575") %>%
  image_annotate("-",size = 60, location = "+1905+740") %>%
  image_annotate("-",size = 60, location = "+1905+905") %>%
  image_annotate("-",size = 60, location = "+1905+1070") %>%
  image_annotate("-",size = 60, location = "+1905+1235")

main_figure<-image_append(c(image_append(c(a,b)),image_append(c(c,d))),stack =T)
image_write(main_figure,"F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/main_figure.png")


a1<-image_read("F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/allium_picture/Alineare_gygax_2022_acrop.JPG")
a2<-image_read("F:/InfoFlora/CWR/Analyses/Resultats_200428/results_spRichness/allium_picture/Alineare_gygax_2022_bcrop.JPG")
a2<-image_resize(a2,'x3456')
a<-c(a1,a2) %>% 
  image_scale('x500') %>% 
  image_append
a<-image_append(c(image_blank(121,500,'white'),
                  a,
                  image_blank(121,500,'white'))) %>% 
  image_annotate('a)',location  = "+30+15",size  =30)


alineare<-image_read("F:/InfoFlora/CWR/Analyses/distributions/1002630_Allium_lineare.png")

b<-image_crop(alineare,"890x700+50+65") %>% 
  image_scale('x500')
b<-image_append(
  c(
    image_blank(119,500,'white'),
    b,
    image_blank(100,500,'white')
  )
) %>% 
  image_annotate('b)', location  = "+30+15",size  =30)


c_fig<-image_crop(alineare,"750x500+1040+225")
c_legend<-image_crop(alineare,"25x870+1815+25") %>% 
  image_resize('45x500!')
c_text<-image_blank(60,500,'white') %>% 
  image_annotate('800',location = '+5+65',size = 25) %>% 
  image_annotate('600',location = '+5+170',size = 25) %>% 
  image_annotate('400',location = '+5+275',size = 25) %>% 
  image_annotate('200',location = '+5+376',size = 25) %>% 
  image_annotate('0',location = '+5+470',size = 25) 
c<-image_append(c(c_fig,c_legend,c_text)) %>% 
  image_annotate('c)', location  = "+30+15",size  =30)

d_fig<-image_crop(alineare,"750x500+2025+1180")
d_legend<-image_crop(alineare,"32x60+2790+1010") 
d_text<-image_blank(110,60,'white') %>% 
  image_annotate('absence',location ="+0+5",size = 20) %>% 
  image_annotate('presence',location ="+0+25",size = 20)
d_legend<-image_append(c(d_legend,d_text)) %>% 
  image_resize('160')
d<-image_append(
  c(
    d_fig,
    image_blank(105,500,'white')
  )
) %>% 
  image_composite(d_legend,offset = '+650+50') %>% 
  image_annotate('d)',location = "+30+15",size =30) 


fig2<-image_append(c(image_append(c(a,b)),image_append(c(c,d))),stack =T)
image_write(fig2,path = "F:/InfoFlora/CWR/Analyses/Resultat_220906/fig2.png")

