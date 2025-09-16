# CU25 BDNF-axis stats (R)
library(readr); library(dplyr); library(ggplot2); library(tibble); library(stringr)
# scores <- score_bdnf_axes(expr, trkb_df, p75_df, boost_by_evidence = TRUE)
meta <- read_csv("CU25/CU25_sample_info.csv", show_col_types = FALSE)
pick_sample_col <- function(df, sample_ids) { ov <- sapply(names(df), function(c) length(intersect(as.character(df[[c]]), sample_ids))); names(ov)[which.max(ov)] }
pick_group_col <- function(df) { cand <- names(df)[str_detect(names(df), regex("group|condition|treat|arm|cohort", ignore_case = TRUE))]; cand <- c(cand, setdiff(names(df), cand)); for (c in cand) { nun <- dplyr::n_distinct(df[[c]]); if (nun >= 2 && nun <= 6) return(c) }; NULL }
stopifnot(exists("scores"))
sample_col <- pick_sample_col(meta, unique(scores$sample))
group_col <- pick_group_col(meta)
meta2 <- meta %>% transmute(sample = .data[[sample_col]], GroupRaw = if (!is.null(group_col)) .data[[group_col]] else NA_character_) %>%
  mutate(Group2 = case_when(str_detect(GroupRaw, regex("id8", ignore_case = TRUE)) ~ "ID8", str_detect(GroupRaw, regex("ctrl|control|vehicle|veh", ignore_case = TRUE)) ~ "Control", TRUE ~ GroupRaw))
out <- scores %>% left_join(meta2, by = "sample") %>% filter(Group2 %in% c("Control","ID8")) %>% mutate(Group2 = factor(Group2, levels = c("Control","ID8")))
welch_t <- function(x, y) { x<-as.numeric(x); y<-as.numeric(y); nx<-length(x); ny<-length(y); mx<-mean(x); my<-mean(y); vx<-var(x); vy<-var(y); se<-sqrt(vx/nx+vy/ny); t<-(mx-my)/se; df<-(vx/nx+vy/ny)^2/((vx^2)/((nx^2)*(nx-1))+(vy^2)/((ny^2)*(ny-1))); p<-2*pt(abs(t), df=df, lower.tail=FALSE); ci<-(mx-my)+qt(c(0.025,0.975),df)*se; s_pooled<-sqrt(((nx-1)*vx+(ny-1)*vy)/(nx+ny-2)); d<-(mx-my)/s_pooled; J<-1-(3/(4*(nx+ny)-9)); g<-d*J; c(mean_diff=mx-my,t=t,df=df,p=p,ci_low=ci[1],ci_high=ci[2],hedges_g=g) }
res <- bind_rows(as.data.frame(t(welch_t(out$trkb_score[out$Group2=="ID8"], out$trkb_score[out$Group2=="Control"]))) %>% mutate(metric="trkb_score"),
                 as.data.frame(t(welch_t(out$p75_score[out$Group2=="ID8"], out$p75_score[out$Group2=="Control"]))) %>% mutate(metric="p75_score"),
                 as.data.frame(t(welch_t(out$bdnf_bias[out$Group2=="ID8"], out$bdnf_bias[out$Group2=="Control"]))) %>% mutate(metric="bdnf_bias")) %>%
  mutate(n_Control=sum(out$Group2=="Control"), n_ID8=sum(out$Group2=="ID8")) %>%
  select(metric, n_Control, n_ID8, mean_diff, t, df, p, ci_low, ci_high, hedges_g) %>% mutate(q_BH=p.adjust(p, method="BH"))
readr::write_csv(res, "bdnf/CU25_bdnf_axis_welch_tests.csv")
