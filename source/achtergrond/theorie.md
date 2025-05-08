# Steekproefgrootte en kleinste detecteerbare effect

## Concept

Het kleinste detecteerbare effect is het kleinste effect waarvoor we, wanneer het zich in werkelijkheid zou voordoen, een bepaalde kans (de power) hebben om het te detecteren.
In dit rapport gaan we uit van een gewenste power van 90% bij een kans op een Type-I fout van 10%.
We aanvaarden dus 10% kans om ten onrechte een effect te veronderstellen dat er niet is.
In de ecologie is het vaak belangrijker om tijdig een effect te detecteren dan om een effect dat er niet is te verwerpen.
Beter op tijd ingrijpen terwijl het (nog) niet nodig is, dan te laat ingrijpen omdat we nog niet helemaal zeker van ons stuk zijn.

Enkel voor eenvoudige experimentele ontwerpen kunnen we de steekproefgrootte analytisch bepalen.
Voor de meeste werkelijke monitoringschema's moeten we terug op simulaties.
Dit houdt in dat we een fictieve dataset genereren volgens het verwachtte ontwerpen en de verwachte eigenschappen.
Vervolgens passen we de statistische analyse toe die we ook op de echte data zouden toepassen.
En dan kijken we of we het effect dat we in de simulatie hebben ingebouwd ook terugvinden.
Uiteraard hangt het effect sterk van de gegeneerde dataset af.
Daarom herhalen we de simulatie een groot aantal keer.
De geschatte power is dan het percentage keren dat we het effect terugvinden.

Zelden hebben we dadelijk de geschikte combinatie van parameters.
Bij het ontwerp van een meetnet gaan we op zoek naar de steekproefgrootte die ons de gewenste power geeft voor het kleinste detecteerbare effect dat we willen vinden.
In deze oefening doen we het omgekeerde: wat is het kleinste detecteerbare effect dat we kunnen vinden met het huidige meetnet?
Het Cmon meetnet is oorspronkelijk uitgewerkt om op basis van 20 jaar data uitspraken te doen over de gemiddelde koolstofvoorraad over de landgebruiken heen.
Momenteel is de vraag welke uitspraken we kunnen doen met Cmon voor de koolstofvoorraad in de bossen voor kortere periodes.
Om dat we voor de natuurherstelverordening elke 6 jaar moeten rapporten, gaan we op zoek naar de kleinste detecteerbare trend over een aantal veelvouden van 6 jaar.
We kijken tevens naar de trend over de volledige periode van 20 jaar.

## Vuistregels

Eens we de steekproefgrootte en het bijhorende kleinste detecteerbare effect voor een bepaalde combinatie van kenmerken kennen, kunnen we aan de hand van vuistregels een redelijke schatting maken voor een andere combinatie.
Deze vuistregels zijn gebaseerd op jaarlijks herhaalde metingen.
Binnen Cmon is gekozen voor een roterend schema waarbij de locaties elke tien jaar opnieuw gemeten worden.
Dat wil zeggen dat we vuistregels met de nodige voorzichtigheid moeten gebruiken wanneer we deze toepassen op een looptijd die korter dan twee keer 10 jaar is.

### Het kleinste detecteerbare effect is rechtevenredig met de variantie

Wanneer we er in slagen om de (totale) variantie te halveren, dan halveert het kleinste detecteerbare effect.

### Het kleinste detecteerbare effect is omgekeerd evenredig met het kwadraat van de steekproefgrootte

Wanneer we de steekproefgrootte verviervoudigen, dan halveert het kleinste detecteerbare effect.

### De looptijd is omgekeerd evenredig met de derde macht van de jaarlijkse steekproefgrootte

Wanneer we de looptijd willen halveren en de detecteerbare jaarlijkse wijziging constant willen houden, dan moeten we de jaarlijkse steekproefgrootte vermenigvuldigen met acht.
Bij de analyse is het totaal aantal metingen belangrijk, niet het aantal metingen per jaar.
Halveren we de looptijd, dan moeten we de jaarlijkse steekproefgrootte verdubbelen om hetzelfde totaal aantal metingen te hebben binnen de looptijd.
Daarnaast is de cumulatieve wijziging over de looptijd ook belangrijk.
Na de helft van de periode is de cumulatieve wijziging dan ook nog maar half zo groot als over de volledige periode.
De vorige vuistregel heeft aan dat we dat moeten compenseren door de jaarlijkse steekproefgrootte te verviervoudigen.
En aangezien we de jaarlijkse steekproefgrootte als met twee moesten vermenigvuldigen om het totaal op peil te houden, moeten we de jaarlijkse steekproefgrootte met acht vermenigvuldigen.

## Analyse van de Cmon gegevens per landgebruik

Voor deze concrete oefening beperken we ons tot een analyse van de Cmon gegevens voor elk landgebruik afzonderlijk.
$OC_{it}$ is de voorraad organische koolstof in de bodem op locatie $i$ op tijdstip $t$.
We veronderstellen dat de voorraad organische koolstof in de bodem log-normaal verdeeld is met gemiddelde $\mu_{it}$ en meetvariantie $\sigma^2_e$ (@eq-lognormal).
De voorraad in alle locaties neemt volgens dezelfde trend $\beta_t$ toe of af (@eq-trend).
De voorraad is locatie $i$ bedraagt gemiddeld een factor $\exp{b_i}$ meer of minder dan de gemiddelde voorraad.
De spreiding tussen de locaties hangt af van de landgebruikvariantie $\sigma^2_l$ (@eq-random).

$$\log(OC_{it}) \sim \mathcal{N}(\mu_{it}, \sigma^2_e)$$ {#eq-lognormal}

$$\mu_{it} = \mu_{i0} + \beta t + b_i$$ {#eq-trend}

$$b_i \sim \mathcal{N}(0, \sigma^2_l)$$ {#eq-random}

We schatten de parameters van het model met behulp van het `glmmTMB` package [@glmmTMB] in de statistische software R [@R].
