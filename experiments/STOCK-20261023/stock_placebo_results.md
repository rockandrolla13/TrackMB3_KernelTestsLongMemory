# Stock trade-duration placebo

Exploratory placebo; process independence is unverified; days are reused.

| Pairing | n | Valid pairs | HSIC wild | HSIC shift | HSIC sub | dCov wild | dCov sub |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| AAPL-MSFT | 500 | 500 | 9.4% | 4.4% | 10.0% | 9.0% | 9.8% |
| AAPL-MSFT | 1000 | 500 | 7.0% | 3.8% | 9.4% | 7.6% | 9.0% |
| AAPL-MSFT | 2000 | 500 | 8.6% | 5.6% | 8.4% | 8.4% | 8.8% |
| AAPL-NVDA | 500 | 500 | 11.0% | 6.8% | 11.6% | 8.2% | 10.2% |
| AAPL-NVDA | 1000 | 500 | 8.4% | 6.0% | 10.0% | 9.0% | 10.0% |
| AAPL-NVDA | 2000 | 500 | 5.8% | 4.8% | 9.2% | 6.0% | 7.2% |
| MSFT-NVDA | 500 | 500 | 9.2% | 5.4% | 10.8% | 9.2% | 11.8% |
| MSFT-NVDA | 1000 | 500 | 7.8% | 5.6% | 10.6% | 9.4% | 10.0% |
| MSFT-NVDA | 2000 | 500 | 7.2% | 5.4% | 10.2% | 6.4% | 11.2% |
| AAPL-AAPL | 500 | 500 | 8.6% | 6.0% | 9.4% | 8.2% | 8.6% |
| AAPL-AAPL | 1000 | 500 | 8.8% | 6.8% | 10.0% | 8.0% | 9.0% |
| AAPL-AAPL | 2000 | 500 | 4.8% | 3.4% | 7.2% | 5.0% | 7.6% |

Each window stays within one session, and paired days are at least 20 trading sessions apart.
Reused days and unknown process independence prevent false-alarm interpretation.
