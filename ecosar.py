import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Dados de exemplo
df = pd.read_excel('Ecosar_data.xlsx')

# Criação do gráfico com tamanho ajustado
fig, ax = plt.subplots(figsize=(12, 8))


# Adicionando faixas coloridas ao gráfico
ax.axhspan(-2.0, 0.0, facecolor='#F4AC6C', alpha=0.5)
ax.axhspan(0.0, 1.0, facecolor='#F7E1BF', alpha=0.5)
ax.axhspan(1.0, 2.0, facecolor='#ADD8E6', alpha=0.5)
ax.axhspan(2.0, 8.0, facecolor='#F08080', alpha=0.5)

# Adicionando texto ao gráfico
ax.text(0, -1, 'Very toxic')
ax.text(0, 0.5, 'Toxic')
ax.text(0, 1.5, 'Harmful')
ax.text(0, 5, 'Not harmful')


ax.scatter(df['Moleculas'], np.log10(df['Fish']), color='blue', marker='o', label='Fish')
ax.scatter(df['Moleculas'], np.log10(df['Daphnid']), color='red', marker='s', label='Daphnid')
ax.scatter(df['Moleculas'], np.log10(df['Green Algae']), color='green', marker='^', label='Green Algae')


# Adicionando a legenda
plt.legend()

# Adicionando a legenda do eixo y com um número subscrito
plt.ylabel(r'log IC$_{50}$ or EC$_{50}$ (mg/L)')
plt.ylim(-2.0, 8.0)
plt.yticks(np.arange(-2.0, 9.0, 1.0))
#ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Exibindo o gráfico
plt.show()
