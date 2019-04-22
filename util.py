import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import pingouin as pg
import rpy2.robjects as robjects
import seaborn as sns
import warnings

from rpy2.rinterface import RRuntimeWarning
from rpy2.robjects import pandas2ri

import effect_size as effect

# set ggplot style
plt.style.use('ggplot')
mpl.rcParams.update({'font.size': 24})

sns.set_style("whitegrid")

# Permite manipular o R em "alto nível"

# Filtro os warning do R


warnings.filterwarnings("ignore", category=RRuntimeWarning)

# Usar pandas

pandas2ri.activate()


def require_r_package(package_name):
    devtools = "devtools::install_github(\"b0rxa/scmamp\")"
    pkg = f"install.packages(\"{package_name}\", repos='https://cran-r.c3sl.ufpr.br/')"
    msg = f"\"You are missing the package '{package_name}', we will now try to install it...\""
    a = """if(!require({0})){{
       print({1})
       {2}
       require({0})
    }}""".format(package_name, msg, devtools if package_name == "scmamp" else pkg)

    return a


class kruskal_test(object):
    def __init__(self, df, val_col: str, group_col: str, sort=True):
        self.df = df
        self.val_col = val_col
        self.group_col = group_col

        if not sort:
            self.df = df.sort_values(by=group_col)

        self.r_dataframe = pandas2ri.py2ri(self.df)

    def apply(self, alpha=0.05, plot=True, filename="kruskal", use_latex=False):
        kruskal = pg.kruskal(dv=self.val_col, between=self.group_col, data=self.df)
        pvalue = kruskal['p-unc'][0]

        if plot:
            chi_squared = kruskal['H'][0]
            degree_freed = kruskal['ddof1'][0]

            p = "< 0.001" if pvalue < 0.001 else (
                "< 0.01" if pvalue < 0.01 else ("< 0.05" if pvalue < 0.05 else (round(pvalue, 3))))

            plt.figure(figsize=(70, 8))
            sns.boxplot(x=self.group_col, y=self.val_col, data=self.df)
            # Jittered BoxPlots
            sns.stripplot(x=self.group_col, y=self.val_col, data=self.df, size=4, jitter=True, edgecolor="gray")
            # Add mean and median lines
            plt.axhline(y=self.df[self.val_col].mean(), color='r', linestyle='--', linewidth=1.5)
            plt.axhline(y=self.df[self.val_col].median(), color='b', linestyle='--', linewidth=2)

            plt.title("")
            plt.suptitle("")
            plt.xlabel(f"\nKruskal-Wallis chi-squared = {chi_squared}, df = {degree_freed}, p = {p}", labelpad=20)
            plt.ylabel('')
            plt.savefig(filename + ('.pgf' if use_latex else '.pdf'), bbox_inches='tight')
            plt.clf()

        # If the Kruskal-Wallis test is significant, a post-hoc analysis can be performed
        # to determine which levels of the independent variable differ from each other level.
        if pvalue < alpha:            
            eff = effect_size(self.df, self.val_col, self.group_col)
            return kruskal, [self._post_hoc_nemenyi(), self._kruskal_multiple_comparisons(), eff.VD_A()]

        return kruskal, None

    def _post_hoc_nemenyi(self):
        """
        Nemenyi test for multiple comparisons
        Zar (2010) suggests that the Nemenyi test is not appropriate for groups with unequal numbers of observations.
        :return: 
        """

        rposthoc = """kruskal.post.hoc <- function(value, group, alpha = 0.05){{
                {0}
                {1}
                {2}                            
                group <- as.factor(group)
                post.results <- tryCatch({{ 
                  kwAllPairsNemenyiTest(value, group, "Tukey")
                }}, warning=function(w) {{      
                  kwAllPairsNemenyiTest(value, group, "Chisquare")
                }})
                
                print("============================= POST-HOC TESTS =============================")                
                
                ## LaTeX formated: Significances highlighted in bold
                ## Pay attention! Read from left to right the comparations. See the TRUE values!!!
                no.diff <- post.results$p.value < alpha
                no.diff[is.na(no.diff)] <- FALSE
                writeTabular(table=post.results$p.value, format='f', bold=no.diff,hrule=0,vrule=0)
                                
                return(post.results)
            }}
            """.format(require_r_package("PMCMRplus"), require_r_package("devtools"), require_r_package("scmamp"))

        kruskal_post_hoc_nemenyi = robjects.r(rposthoc)

        return kruskal_post_hoc_nemenyi(self.r_dataframe.rx2(self.val_col), self.r_dataframe.rx2(self.group_col))

    def _kruskal_multiple_comparisons(self):
        """
        Multiple comparison test between treatments or treatments versus control after Kruskal-Wallis test.

        This test helps determining which groups are different with pairwise comparisons adjusted
        appropriately for multiple comparisons. Those pairs of groups which have observed differences
        higher than a critical value are considered statistically different at a given significance level.
        :return:
        """

        rmulti = """kruskal.mc <- function(value, group, data, alpha = 0.05){{
            {0}
            group <- as.factor(group)
            return(kruskalmc(value ~ group, probs=alpha, data=data, cont=NULL))
        }}
        """.format(require_r_package("pgirmess"))

        kruskal_mc = robjects.r(rmulti)

        return kruskal_mc(self.r_dataframe.rx2(self.val_col), self.r_dataframe.rx2(self.group_col), self.r_dataframe)


class effect_size(object):
    def __init__(self, df, val_col: str, group_col: str, sort=True):
        self.df = df
        self.val_col = val_col
        self.group_col = group_col

        if not sort:
            self.df = df.sort_values(by=group_col)

        print("EFFECT")
        XX = df.groupby([group_col]).agg(['count'])
        print(XX)

    def VD_A(self):
        return effect.VD_A_DF(self.df, self.val_col, self.group_col)


class friedman_test(object):
    def __init__(self, filename, val_col: str, group_col: str, instance_col: str):
        """

        :param df: pandas DataFrame object
        :param val_col: Name of the column that contains values (fitness values).
        :param group_col: Name of the column that contains group names (algorith namess).
        :param instance_col: Name of the column that contains instance names (problem names).
        """

        self.df = pd.read_csv(filename, ";")
        self.filename = filename
        self.val_col = val_col
        self.group_col = group_col
        self.instance_col = instance_col

        self.r_dataframe = pandas2ri.py2ri(self.df)

    def apply(self):
        return pg.friedman(dv=self.val_col, within=self.group_col, subject=self.instance_col, data=df)

    def apply_r(self, instance_col: str = "Instance"):
        rfriedman = """friedman.post.hoc <- function(filename, instance_col, alpha = 0.05){{
            {0}
            {1}
            {2}
            #read the data
            data <- read.csv(filename, sep=";")
            #remove the Instance column
            data <- data[,!(names(data) %in% c(instance_col))]
            data <- data.matrix(data)
            
            pre.results <- friedmanTest(data)
            
            if(pre.results$p.value < alpha){{
                post.results <- posthoc.friedman.nemenyi.test(data)
                friedman.mc <- friedmanmc(data)
                setEPS()
                postscript("{3}_plotCD.eps", width = 6.0, height = 3.0)
                output <- scmamp::plotCD(data, alpha=0.05 )
                dev.off()
                
                list.to.return <- list(Friedman = pre.results, PostHoc = post.results, FriedmanMC = friedman.mc)
                return(list.to.return)                
            }}else{{
                list.to.return <- list(Friedman = pre.results, PostHoc = NULL, FriedmanMC = NULL)
                return(list.to.return)
            }}
        }}    
        """.format(require_r_package("PMCMRplus"), require_r_package("PMCMR"), require_r_package("pgirmess"),
                   self.filename.replace(".csv", ""))

        friedman_r = robjects.r(rfriedman)

        return friedman_r(self.filename, instance_col)


if __name__ == '__main__':
    df = pd.read_csv("results/deeplearning4j@deeplearning4j/stats_tcfail.csv", ";")
    k = kruskal_test(df, 'fitness', 'policy')
    kruskal, posthoc = k.apply()

    if posthoc is not None:
        print(posthoc[0])
        print(posthoc[1])
        print(posthoc[2])

    f = friedman_test("results/Friedman_Test.csv", 'fitness', 'algorithm', 'instance')
    print(f.apply_r())
