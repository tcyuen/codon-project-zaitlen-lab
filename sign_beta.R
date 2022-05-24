##sign of beta

add_sign <- read_tsv("height_data_percent_change.tsv", col_names = TRUE)

#change value of beta 
add_sign$beta <- ifelse(add_sign$ref != add_sign$ancestral_allele,
                                    -add_sign$beta, add_sign$beta)

#add column to see if beta is positive or negative
sign_of_beta <- as.factor(ifelse(add_sign$beta < 0, "-", "+"))

add_sign$beta_sign <- sign_of_beta

#sign_model <- lm(beta_sign ~ delta, data = add_sign)

sign_model <- lm(delta ~ beta_sign, data = add_sign)

summary(sign_model)$coef

contrasts(add_sign$beta_sign)

logit <- glm(beta_sign ~ delta, data = add_sign, family = "binomial")
summary(logit)


## run on p-value thresholds
add_sign$pval <- as.numeric(as.character(add_sign$pval))

beta_signs <- add_sign[add_sign$pval < 0.05, ]$beta_sign

delta_codon_freqs <- add_sign[add_sign$pval < 0.05, ]$delta

summary(glm(beta_signs ~ delta_codon_freqs, family = "binomial"))



#pval <0.01
beta_signs <- add_sign[add_sign$pval < 0.01, ]$beta_sign

delta_codon_freqs <- add_sign[add_sign$pval < 0.01, ]$delta


summary(glm(beta_signs ~ delta_codon_freqs, family = "binomial"))


#pval <0.001
beta_signs <- add_sign[add_sign$pval < 0.001, ]$beta_sign

delta_codon_freqs <- add_sign[add_sign$pval < 0.001, ]$delta


summary(glm(beta_signs ~ delta_codon_freqs, family = "binomial"))

#pval < 0.0001
beta_signs <- add_sign[add_sign$pval < 0.0001, ]$beta_sign

delta_codon_freqs <- add_sign[add_sign$pval < 0.0001, ]$delta

summary(glm(beta_signs ~ delta_codon_freqs, family = "binomial"))


#pval < 0.00001
beta_signs <- add_sign[add_sign$pval < 0.00001, ]$beta_sign

delta_codon_freqs <- add_sign[add_sign$pval < 0.00001, ]$delta

summary(lm(beta_signs ~ delta_codon_freqs))
sign_logit <- glm(beta_signs ~ delta_codon_freqs)

ggplotRegression(sign_logit)

#pval < 0.000001
beta_signs <- add_sign[add_sign$pval < 0.000001, ]$beta_sign

delta_codon_freqs <- add_sign[add_sign$pval < 0.000001, ]$delta

summary(lm(beta_signs ~ delta_codon_freqs))
sign_logit <- glm(beta_signs ~ delta_codon_freqs)

ggplotRegression(sign_logit)

#pval < 0.0000001
beta_signs <- add_sign[add_sign$pval < 0.0000001, ]$beta_sign

delta_codon_freqs <- add_sign[add_sign$pval < 0.0000001, ]$delta

summary(lm(beta_signs ~ delta_codon_freqs))
sign_logit <- glm(beta_signs ~ delta_codon_freqs)

ggplotRegression(sign_logit)

