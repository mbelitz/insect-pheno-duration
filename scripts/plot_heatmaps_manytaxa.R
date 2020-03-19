library(rgbif)
library(dplyr)
library(lubridate)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(ridigbio)

heat_map_obs <- function(binomial){
  
  idig <- idig_search_records(limit = 10000, rq = list(scientificname = binomial), 
                              fields = c("data.dwc:decimalLongitude", "data.dwc:decimalLatitude",
                                         "data.dwc:year", "data.dwc:month", "data.dwc:eventDate",
                                         "data.dwc:scientificName", "data.dwc:basisOfRecord")) 
  
  names(idig) <- c("decimalLongitude", "decimalLatitude", "year", "month", "eventDate",
                   "scientificName", "basisOfRecord")
  
  #rgbif
  gbif <- occ_search(scientificName = binomial,
                     hasCoordinate = TRUE, country = "US",
                     limit = 200000)
  
  gbif_df <- gbif$data %>% 
    filter(basisOfRecord == "HUMAN_OBSERVATION") %>% 
    select(decimalLongitude, decimalLatitude, year, month, eventDate, 
           scientificName, basisOfRecord)
  
  total <- rbind(gbif_df, idig) %>% 
    filter(!is.na(year)) %>% 
    filter(!is.na(decimalLatitude)) %>% 
    filter(!is.na(decimalLongitude)) %>% 
    filter(!is.na(eventDate))
  
  total$decimalLongitude <- as.numeric(total$decimalLongitude)
  total$decimalLatitude <- as.numeric(total$decimalLatitude)
  total$year <- as.numeric(total$year)
  
  total <- total %>% 
    mutate(decade = floor(year/10)*10)
  
  us <- rnaturalearth::ne_countries(country = "United States of America", 
                                    returnclass = "sf")
  
 a <-  ggplot() + 
    geom_sf(us, mapping = aes()) +
    geom_bin2d(total,mapping = aes(x = decimalLongitude, y = decimalLatitude), alpha = 0.9, binwidth = c(2,2)) + 
    coord_sf( xlim = c(-125, -65), ylim = c(20, 60)) +
    scale_fill_continuous(type = "viridis",trans = "log10") + 
    ggtitle(binomial) + 
    facet_wrap(~decade)

 ggsave(filename = paste("map_outputs/",binomial,".jpeg",sep = ""), dpi = 100, plot = a)
 
}

saturniidae <- list('Automeris io',
              'Eacles imperialis',
              'Antheraea polyphemus',
              'Dryocampa rubicunda',
              'Hemileuca eglanterina',
              'Actias luna',
              'Callosamia promethea',
              'Hyalophora cecropia',
              'Hemileuca maia',
              'Anisota stigma')

lapply(X = moths,FUN = heat_map_obs)

cicadas <- list(
'Platypedia putnami',
'Magicicada septendecim',
'Neotibicen tibicen',
'Megatibicen auletes',
'Pacarina puella',
'Diceroprocta apache',
'Megatibicen pronotalis',
'Diceroprocta eugraphica',
'Okanagana bella')

lapply(cicadas, heat_map_obs)

mayflies <- list('Stenacron interpunctatum',
'Hexagenia limbata',
'Caenis latipennis',
'Stenonema femoratum',
'Baetis tricaudatus',
'Maccaffertium terminatum',
'Ephemerella invaria',
'Maccaffertium vicarium',
'Baetis flavistriga',
'Hexagenia bilineata')

lapply(mayflies, heat_map_obs)

dragonflies <- list('Ischnura verticalis',
                    'Calopteryx maculata',
                    'Pachydiplax longipennis',
                    'Enallagma boreale',
                    'Hetaerina americana',
                    'Lestes disjunctus',
                    'Enallagma civile',
                    'Sympetrum obtrusum',
                    'Argia fumipennis',
                    'Erythemis simplicicollis')


lapply(X = dragonflies, FUN = heat_map_obs)
heat_map_obs("Erythemis simplicicollis")

dobsonflies <- list('Nigronia serricornis',
                    'Corydalus cornutus',
                    'Chauliodes pectinicornis',
                    'Corydalus luteus',
                    'Chauliodes rastricornis',
                    'Nigronia fasciata',
                    'Corydalus peruvianus',
                    'Platyneuromus soror',
                    'Sialis americana',
                    'Neohermes concolor'
)

lapply(dobsonflies, heat_map_obs)

scarabs <- list('Tomarus gibbosus',
'Popillia japonica',
'Onthophagus hecate',
'Canthon imitator',
'Onthophagus gazella',
'Oxygrylius ruginasus',
'Ateuchus histeroides',
'Onthophagus pennsylvanicus',
'Canthon indigaceus',
'Ataenius platensis')

lapply(scarabs, heat_map_obs)

silphidae <- list('Nicrophorus orbicollis',
                  'Nicrophorus tomentosus',
                  'Necrophila americana',
                  'Nicrophorus guttula',
                  'Thanatophilus lapponicus',
                  'Nicrophorus carolinus',
                  'Nicrophorus mexicanus',
                  'Nicrophorus marginatus',
                  'Oiceoptoma noveboracense',
                  'Nicrophorus defodiens'
)

lapply(silphidae, heat_map_obs)

ladybugs <- list(
'Hippodamia convergens',
'Coccinella transversoguttata',
'Coccinella novemnotata',
'Coleomegilla maculata',
'Harmonia axyridis',
'Adalia bipunctata',
'Psyllobora vigintimaculata',
'Coccinella trifasciata',
'Hippodamia quinquesignata',
'Hippodamia parenthesis')

lapply(ladybugs, heat_map_obs)

wasps <- list(
'Bembix americana',
'Philanthus gibbosus',
'Plenoculus davisi',
'Trypoxylon lactitarse',
'Oxybelus emarginatus',
'Microbembex monodonta',
'Sphecius speciosus',
'Bembix amoena')

lapply(wasps, heat_map_obs)

spiderwasps <- list(
'Episyron vagabundus',
'Paracyphononyx funereus',
'Priocnessus nebulosus',
'Poecilopompilus familiaris',
'Sericopompilus apicalis',
'Aporus niger',
'Dipogon sericeus',
'Caliadurgus fasciatellus',
'Ceropales maculata',
'Minagenia osoria')

lapply(spiderwasps, heat_map_obs)

parasiticwasps <- list(
'Diplazon laetatorius',
'Trogomorpha trogiformis',
'Barichneumon neosorex',
'Itoplectis conquisitor',
'Pimpla pedalis',
'Zaglyptus pictilis')

lapply(parasiticwasps, heat_map_obs)

sphinx <- list(
'Smerinthus cerisyi',
'Hemaris diffinis',
'Paonias excaecata',
'Paonias myops',
'Smerinthus jamaicensis',
'Hyles lineata',
'Darapsa myron',
'Hemaris thysbe',
'Sphinx vashti',
'Manduca sexta')

lapply(sphinx, heat_map_obs)

carabidae <- list(
  'Cicindela oregona',
  'Cicindela punctulata',
  'Pterostichus adstrictus',
  'Calathus ruficollis',
  'Cicindela repanda',
  'Cicindela tranquebarica',
  'Dicheirus piceus',
  'Pterostichus lama',
  'Cicindela scutellaris'
)

lapply(carabidae, heat_map_obs)

stagbeetles <- list(
  'Lucanus capreolus',
  'Lucanus placidus',
  'Sinodendron rugosum',
  'Ceruchus piceus',
  'Lucanus mazama',
  'Platycerus virescens',
  'Platycerus oregonensis',
  'Pseudolucanus capreolus',
  'Dorcus parallelus',
  'Ceruchus punctatus'
)

lapply(stagbeetles, heat_map_obs)

longhorn <- list(
'Gnathacmaeops pratensis',
'Xylotrechus colonus',
'Xylotrechus sagittatus',
'Neoclytus acuminatus',
'Prionus californicus',
'Megacyllene robiniae',
'Acmaeops proteus',
'Clytus ruricola',
'Mallodon dasystomum')

lapply(longhorn, heat_map_obs)

metallicbeet <- list(
  'Acmaeodera rubronotata',
  'Acmaeodera mixta',
  'Acmaeodera gibbula',
  'Acmaeodera amplicollis',
  'Acmaeodera decipiens',
  'Acmaeodera bowditchi',
  'Brachys aerosus',
  'Acmaeodera pulchella',
  'Acmaeodera scalaris',
  'Anthaxia inornata'
)

lapply(metallicbeet, heat_map_obs)

robberflies <- list(
  'Holcocephala abdominalis',
  'Leptogaster flavipes',
  'Psilonyx annulatus',
  'Holcocephala fusca',
  'Holcocephala calva',
  'Tipulogaster glabrata',
  'Lasiopogon terricolus',
  'Leptogaster obscuripennis',
  'Megaphorus pulchrus'
)

lapply(robberflies, heat_map_obs)

horseflies <- list(
  'Hybomitra lasiophthalma',
  'Tabanus quinquevittatus',
  'Chrysops mitis',
  'Hybomitra epistates',
  'Chrysops vittatus',
  'Tabanus abactor',
  'Chrysops univittatus',
  'Hybomitra illota',
  'Chrysops excitans',
  'Chrysops indus'
)

lapply(horseflies, heat_map_obs)

clickbeetles <- list(
  'Limonius aeger',
  'Horistonotus simplex',
  'Agrypnus rectangularis',
  'Limonius auripilis',
  'Hypnoidus bicolor',
  'Melanotus similis',
  'Conoderus bellus',
  'Limonius basilaris',
  'Conoderus aversus',
  'Melanotus morosus'
)

lapply(clickbeetles, heat_map_obs)

blisterbeet <- list(
'Epicauta bispinosa',
'Epicauta ferruginea',
'Epicauta pensylvanica',
'Epicauta maculata',
'Epicauta pardalis',
'Lytta cyanipennis',
'Epicauta puncticollis',
'Nemognatha lutea',
'Epicauta polingi',
'Epicauta normalis')

lapply(blisterbeet, heat_map_obs)

water <- list(
'Cercyon quisquilius',
'Tropisternus collaris',
'Enochrus ochraceus',
'Tropisternus affinis',
'Tropisternus lateralis',
'Cercyon fimbriatus',
'Enochrus pygmaeus',
'Enochrus hamiltoni',
'Berosus exiguus',
'Crenitis snoqualmie')

lapply(water, heat_map_obs)

leafbeet <- list(
  'Diabrotica undecimpunctata',
  'Gastrophysa cyanea',
  'Diabrotica adelpha',
  'Gratiana pallidula',
  'Leptinotarsa decemlineata',
  'Chrysomela scripta',
  'Diachus auratus',
  'Tricholochmaea cavicollis',
  'Chrysomela interrupta'
)

lapply(leafbeet, heat_map_obs)


bees <- list(
  'Bombus impatiens',
  'Bombus bifarius',
  'Bombus flavifrons',
  'Bombus pensylvanicus',
  'Bombus occidentalis',
  'Bombus mixtus',
  'Bombus sylvicola',
  'Bombus vagans')

lapply(bees, heat_map_obs)

heat_map_obs("Apis mellifera")

